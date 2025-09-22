program simple_network_example_p
    use hyperSIS_kinds_mod, only: dp, i4, fmt_general
    use hyperSIS_network_mod, only : network_t
    use hyperSIS_dynamics_mod, only: dyn_parameters_t, net_state_base_t, net_state_choose
    use hyperSIS_program_common_mod, only: read_network, set_dyn_params, check_qs_method, proc_net_state_gen, proc_export_states, check_export_nodes_and_edges_state, set_initial_number_of_infected_nodes

    use datastructs_mod, only: measure_controller_t, statistical_measure_t

    use rndgen_mod, only: rndgen

    use datastructs_mod, only: log_write, set_level, set_verbose, LOGGER_OK, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG
    implicit none

    ! Network and dynamics variables
    type(network_t) :: net
    type(dyn_parameters_t) :: dyn_params
    class(net_state_base_t), allocatable :: dynamics_state
    procedure(proc_net_state_gen), pointer :: after_dynamics_step => null()
    procedure(proc_export_states), pointer :: export_states => null()

    ! Random number generator
    type(rndgen) :: dyn_gen

    ! Auxiliary variables
    integer(kind=i4) :: i, i_sample

    ! Parameters (via CLI)
    logical :: logger_verbose
    character(len=:), allocatable :: logger_level

    integer(kind=i4) :: rnd_seed
    integer(kind=i4) :: n_samples
    integer(kind=i4) :: initial_number
    character(len=:), allocatable :: output_prefix, edges_file, algorithm, sampler_choice, time_scale
    real(kind=dp) :: par_b, par_theta, initial_infected_fraction, beta_1, tmax
    logical :: use_qs, use_example_network, export_states_flag

    ! Measurers
    type(measure_controller_t) :: time_control
    type(statistical_measure_t) :: time_average
    type(statistical_measure_t) :: rho_average

    ! Get input data
    call handle_cli()

    ! Initialize main structures (network, parameters, measurers)
    call init_main()

    ! Loop over dynamics samples
    do i_sample = 1, n_samples
        ! Allocate dynamics state
        call init_dynamics(i_sample + (n_samples + 1))

        ! Run the dynamics up to time tmax, while there are infected nodes
        call loop_over_time()

        ! Write results to a file
        call write_results()

        ! Deallocate dynamics state for the next sample
        deallocate(dynamics_state)
    end do

contains

    subroutine handle_cli()
        use flap, only : command_line_interface
        type(command_line_interface) :: cli    ! Command Line Interface (CLI)
        character(len=2048) :: cli_string
        integer(kind=i4) :: cli_error

        call cli%init(progname    = 'Run temporal dynamics on a hypergraph network', &
            description = 'This program runs a temporal dynamics on a hypergraph network, using the Gillespie algorithm. The network can be read from a file or generated as a small example. The dynamics can be configured with various parameters.', &
            version     = '1.0', &
            authors      = 'Wesley Cota')

        ! Add options to the CLI
        ! IO parameters
        call cli%add(switch='--output', &
            switch_ab='-o', &
            help='Output prefix for result files', &
            required=.false., &
            act='store', &
            def='./output', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --output'
        call cli%add(switch='--remove-files', &
            switch_ab='-rm', &
            help='Remove existing output files before running', &
            required=.false., &
            act='store', &
            def='false', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --remove-files'
        call cli%add(switch='--edges-file', &
            switch_ab='-e', &
            help='File containing the edge list of the network as input', &
            required=.false., &
            act='store', &
            def='', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --edges-file'

        ! Dynamics and sampler
        call cli%add(switch='--algorithm', &
            switch_ab='-a', &
            help='Dynamics algorithm to use', &
            required=.false., &
            act='store', &
            choices='HB_OGA,NB_OGA', &
            def='HB_OGA', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --algorithm'
        call cli%add(switch='--sampler', &
            switch_ab='-s', &
            help='Sampler choice', &
            required=.false., &
            act='store', &
            choices='rejection_maxheap,btree', &
            def='rejection_maxheap', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --sampler'

        ! Temporal parameters
        call cli%add(switch='--tmax', &
            switch_ab='-t', &
            help='Maximum simulation time', &
            required=.false., &
            act='store', &
            def='1000.0', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --tmax'
        call cli%add(switch='--use-qs', &
            switch_ab='-qs', &
            help='Use Quasi-Stationary method', &
            required=.false., &
            act='store', &
            def='false', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --use-qs'

        ! Sampling and statistics
        call cli%add(switch='--n-samples', &
            switch_ab='-ns', &
            help='Number of samples to average over', &
            required=.false., &
            act='store', &
            def='10', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --n-samples'
        call cli%add(switch='--time-scale', &
            switch_ab='-ts', &
            help='Time scale to use', &
            required=.false., &
            act='store', &
            choices='powerlaw,uniform', &
            def='powerlaw', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --time-scale'

        ! Dynamical parameters
        call cli%add(switch='--initial-fraction', &
            switch_ab='-if', &
            help='Initial fraction of infected nodes', &
            required=.false., &
            act='store', &
            def='1.0', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --initial-fraction'
        call cli%add(switch='--initial-number', &
            switch_ab='-in', &
            help='Initial number of infected nodes (overrides initial-fraction)', &
            required=.false., &
            act='store', &
            def='0', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --initial-number'
        call cli%add(switch='--beta1', &
            switch_ab='-b1', &
            help='Infection rate parameter beta1', &
            required=.true., &
            act='store', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --beta1'
        call cli%add(switch='--par-b', &
            switch_ab='-pb', &
            help='Dynamical parameter b', &
            required=.false., &
            act='store', &
            def='0.5', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --par-b'
        call cli%add(switch='--par-theta', &
            switch_ab='-pt', &
            help='Dynamical parameter theta', &
            required=.false., &
            act='store', &
            def='0.5', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --par-theta'

        ! Additional parameters
        call cli%add(switch='--export-states', &
            switch_ab='-es', &
            help='Export the states of nodes and edges at the end of each sample (be careful with large networks)', &
            required=.false., &
            act='store', &
            def='false', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --export-states'
        call cli%add(switch='--seed', &
            switch_ab='-rs', &
            help='Random seed for the dynamics', &
            required=.false., &
            act='store', &
            def='42', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --seed'
        call cli%add(switch='--verbose', &
            switch_ab='-vv', &
            help='Enable verbose logging', &
            required=.false., &
            act='store', &
            def='true', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --verbose'
        call cli%add(switch='--verbose-level', &
            switch_ab='-vl', &
            help='Logging level', &
            required=.false., &
            act='store', &
            choices='error,warning,info,debug', &
            def='info', &
            error=cli_error); if (cli_error /= 0) error stop 'Error adding --verbose-level'

        call cli%parse(error=cli_error); if (cli_error /= 0) error stop 'Error parsing command line arguments'

        ! Collect parameters from CLI
        call cli%get(switch='--output', val=cli_string, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --output'
        output_prefix = trim(adjustl(cli_string))
        call cli%get(switch='--remove-files', val=export_states_flag, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --remove-files'
        call cli%get(switch='--edges-file', val=cli_string, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --edges-file'
        edges_file = trim(adjustl(cli_string))
        call cli%get(switch='--algorithm', val=cli_string, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --algorithm'
        algorithm = trim(adjustl(cli_string))
        call cli%get(switch='--sampler', val=cli_string, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --sampler'
        sampler_choice = trim(adjustl(cli_string))
        use_example_network = (len_trim(edges_file) == 0)
        call cli%get(switch='--tmax', val=tmax, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --tmax'
        call cli%get(switch='--use-qs', val=use_qs, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --use-qs'
        call cli%get(switch='--n-samples', val=n_samples, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --n-samples'
        call cli%get(switch='--time-scale', val=cli_string, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --time-scale'
        time_scale = trim(adjustl(cli_string))
        call cli%get(switch='--initial-fraction', val=initial_infected_fraction, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --initial-fraction'
        call cli%get(switch='--initial-number', val=initial_number, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --initial-number'
        call cli%get(switch='--beta1', val=beta_1, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --beta1'
        call cli%get(switch='--par-b', val=par_b, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --par-b'
        call cli%get(switch='--par-theta', val=par_theta, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --par-theta'
        call cli%get(switch='--export-states', val=export_states_flag, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --export-states'
        call cli%get(switch='--seed', val=rnd_seed, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --seed'
        call cli%get(switch='--verbose', val=logger_verbose, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --verbose'
        call cli%get(switch='--verbose-level', val=cli_string, error=cli_error); if (cli_error /= 0) error stop 'Error parsing --verbose-level'
        logger_level = trim(adjustl(cli_string))

        call set_verbose(logger_verbose)
        select case (trim(adjustl(logger_level)))
            case ('error')
                call set_level(LOG_ERROR)
            case ('warning')
                call set_level(LOG_WARNING)
            case ('info')
                call set_level(LOG_INFO)
            case ('debug')
                call set_level(LOG_DEBUG)
            case default
                call set_level(LOG_INFO)
        end select

        if (len(edges_file) == 0) then
            use_example_network = .true.
        else
            use_example_network = .false.
        end if

    end subroutine

    subroutine init_main()

        call time_control%init(time_scale)
        call time_average%init(time_control%get_max_array_size(real(tmax,dp)))
        call rho_average%init(time_control%get_max_array_size(real(tmax,dp)))

        select case (use_example_network)
          case (.true.)
            ! Generate a simple network
            call generate_example_network(net)
          case (.false.)
            ! Read the network from a file
            call read_network(net, edges_file)
        end select

        ! Clear and check the network (always necessary)
        call net%clear_and_check_all(min_order=1)

        ! Set initial number of infected nodes if specified
        call set_initial_number_of_infected_nodes(net, initial_infected_fraction, initial_number)

        ! Set the dynamical parameters
        call set_dyn_params(net, dyn_params, par_b, par_theta)

        ! Check if the chosen QS method is compatible with the dynamics
        call check_qs_method(after_dynamics_step, use_qs)

        ! Check if the export states procedure is set
        call check_export_nodes_and_edges_state(export_states, export_states_flag)

    end subroutine init_main

    subroutine init_dynamics(extra_seed)
        integer(kind=i4), intent(in), optional :: extra_seed
        integer(kind=i4) :: seed_offset

        if (present(extra_seed)) then
            seed_offset = extra_seed
        else
            seed_offset = 0
        end if

        ! reset the time control
        call time_control%reset()
        call dyn_gen%init(rnd_seed + seed_offset) ! initialize the random generator with a seed

        ! select dynamics algorithm
        call net_state_choose(dynamics_state, algorithm)

        ! allocate and initializes the dynamics with a random configuration
        call dynamics_state%init_config(net=net, gen=dyn_gen, params=dyn_params, sampler_choice=sampler_choice, fraction=initial_infected_fraction)

        ! Set the beta scale parameter
        ! by default, alpha = 1.0
        dynamics_state%params%beta_scale = beta_1

        ! Initialize the dynamics
        call dynamics_state%dynamics_init(net)
    end subroutine init_dynamics

    subroutine loop_over_time()

        ! main dynamics loop
        do while (dynamics_state%time < tmax)

            ! get Gillespie time
            call dynamics_state%just_update_dt(net, dyn_gen)

            ! update the state
            call dynamics_state%dynamics_step(net, dyn_gen)

            ! do something after the dynamics step
            call after_dynamics_step(net, dynamics_state, dyn_gen)

            ! collect the new data
            call collect_data()

            ! leave if no infected nodes
            if (dynamics_state%get_num_infected() == 0) then
                call log_write(LOG_DEBUG, "No more events can occur, stopping the dynamics at t = ", dynamics_state%time)
                exit ! exit if no infected nodes
            end if
        end do
    end subroutine loop_over_time

    subroutine collect_data()
        integer(kind=i4), allocatable :: time_pos_array(:)
        integer(kind=i4) :: time_pos

        !$omp critical
        time_pos_array = time_control%get_pos_array(dynamics_state%time)

        do time_pos = 1, size(time_pos_array)
            call time_average%add_point(time_pos_array(time_pos), dynamics_state%time)
            call rho_average%add_point(time_pos_array(time_pos), 1.0_dp * dynamics_state%get_num_infected() / net%num_nodes)

            call export_states(net, dynamics_state, get_filename('states_nodes.dat'), get_filename('states_edges.dat'))
        end do
        !$omp end critical

    end subroutine collect_data

    subroutine write_results()
        integer(kind=i4) :: unidade_arquivo
        integer(kind=i4) :: time_pos

        !$omp critical
        ! write time average
        open(newunit=unidade_arquivo, file=get_filename('results.dat'), status='replace', action='write', form='formatted')
        ! write even when not finished
        ! For that, we need to know the maximum number of samples

        write(unidade_arquivo, fmt_general) '# Number of samples: ', time_average%max_n_samples
        write(unidade_arquivo, fmt_general) '# time', 'rho', 'rho_variance', 'n_samples'
        do time_pos = 1, time_control%get_max_array_size(real(tmax,dp))
            if (rho_average%n_samples(time_pos) > 0) then
                ! We use `use_max=.true.` since the missing samples have zero infected nodes
                write(unidade_arquivo, fmt_general) time_average%get_mean(time_pos), rho_average%get_mean(time_pos, use_max=.true.), rho_average%get_variance(time_pos, use_max=.true.), rho_average%n_samples(time_pos)
            end if
        end do
        close(unidade_arquivo)
        !$omp end critical

    end subroutine write_results

    subroutine generate_example_network(net)
        use hyperSIS_network_mod, only : hyperedge, network
        type(network_t), intent(inout) :: net
        integer(kind=i4) :: edge_id, node_pos, node_id

        net = network(num_nodes = 10, num_edges = 7)

        net%edges(1) = hyperedge([1, 2, 3])
        net%edges(2) = hyperedge([2, 3, 4])
        net%edges(3) = hyperedge([4, 5])
        net%edges(4) = hyperedge([1, 5])
        net%edges(5) = hyperedge([6, 7, 8])
        net%edges(6) = hyperedge([1, 2, 3, 7, 8, 9])
        net%edges(7) = hyperedge([9, 10])

        net%nodes(:)%degree = 0

        ! for each edge, collect the degree of each node
        do edge_id = 1, net%num_edges
            do node_pos = 1, net%edges(edge_id)%order + 1
                node_id = net%edges(edge_id)%nodes(node_pos)
                net%nodes(node_id)%degree = net%nodes(node_id)%degree + 1
            end do
        end do

        ! allocate the arrays of edges for each node
        do node_id = 1, net%num_nodes
            allocate(net%nodes(node_id)%edges(net%nodes(node_id)%degree))
            net%nodes(node_id)%degree = 0 ! reset to use as position counter
        end do

        ! for each edge, fill the edges array of each node
        do edge_id = 1, net%num_edges
            do node_pos = 1, net%edges(edge_id)%order + 1
                node_id = net%edges(edge_id)%nodes(node_pos)
                net%nodes(node_id)%degree = net%nodes(node_id)%degree + 1
                net%nodes(node_id)%edges(net%nodes(node_id)%degree) = edge_id
            end do
        end do

    end subroutine generate_example_network

    function get_filename(filename) result(pathname)
        character(len=*), intent(in) :: filename
        character(len=256) :: pathname
        write(pathname, '(g0)') trim(adjustl(output_prefix))//'_'//trim(adjustl(filename))
    end function get_filename

end program simple_network_example_p
