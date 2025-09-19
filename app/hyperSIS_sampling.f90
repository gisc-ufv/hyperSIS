program simple_network_example_p
    use hyperSIS_kinds_mod, only: dp, i4, fmt_general
    use hyperSIS_network_mod, only : network_t
    use hyperSIS_dynamics_mod, only: dyn_parameters_t, net_state_base_t, net_state_choose
    use hyperSIS_program_common_mod, only: read_network, set_dyn_params, check_qs_method, proc_net_state_gen, proc_export_states, check_export_nodes_and_edges_state

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
    integer(kind=i4) :: rnd_seed
    integer(kind=i4) :: n_samples
    character(len=:), allocatable :: edges_file, algorithm, sampler_choice, time_scale, logger_level, output_prefix
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

        call set_verbose(.true.)
        call set_level(LOG_DEBUG)

        ! Default parameters
        use_qs = .true.
        rnd_seed = 9342342
        n_samples = 10
        edges_file = 'example.edgelist'
        par_b = 0.5_dp
        par_theta = 0.5_dp
        initial_infected_fraction = 1.0_dp
        beta_1 = 0.1_dp
        tmax = 1000.0_dp
        time_scale = 'powerlaw' ! or powerlaw
        algorithm = 'HB_OGA' ! or NB-OGA
        sampler_choice = 'rejection_maxheap' ! or btree

        use_example_network = .false. ! if true, generate a small example network; if false, read from edges_file
        export_states_flag = .true. ! if true, export the states of nodes and edges at the end of each sample ! better use with small networks and uniform time scale

        ! Here you can add code to handle command line arguments to override defaults
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

            call export_states(net, dynamics_state, 'nodes.dat', 'edges.dat')
        end do
        !$omp end critical

    end subroutine collect_data

    subroutine write_results()
        integer(kind=i4) :: unidade_arquivo
        integer(kind=i4) :: time_pos

        !$omp critical
        ! write time average
        open(newunit=unidade_arquivo, file='teste', status='replace', action='write', form='formatted')
        ! write even when not finished
        ! For that, we need to know the maximum number of samples

        write(unidade_arquivo, fmt_general) '# Number of samples: ', time_average%max_n_samples
        write(unidade_arquivo, fmt_general) '# time', 'rho', 'rho_variance', 'n_samples'
        do time_pos = 1, time_control%get_max_array_size(real(tmax,dp))
            if (rho_average%n_samples(time_pos) > 0) then
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

end program simple_network_example_p
