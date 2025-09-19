program simple_network_example_p
    use hyperSIS_kinds_mod
    use hyperSIS_network_mod
    use hyperSIS_network_io_mod
    use hyperSIS_dynamics_mod

    use datastructs_mod, only: measure_controller_t, statistical_measure_t

    use rndgen_mod
    implicit none

    ! Network and dynamics variables
    type(network_t) :: net
    type(dyn_parameters_t) :: dyn_params
    class(net_state_base_t), allocatable :: dynamics_state

    ! Random number generator
    type(rndgen) :: dyn_gen

    ! Auxiliary variables
    integer(kind=i4) :: i, i_sample

    ! Parameters
    integer(kind=i4) :: rnd_seed
    integer(kind=i4) :: n_samples
    character(len=:), allocatable :: edges_file, algorithm, sampler_choice, time_scale
    real(kind=dp) :: par_b, par_theta, initial_infected_fraction, beta_1, tmax

    ! Measurers
    type(measure_controller_t) :: time_control
    type(statistical_measure_t) :: time_average
    type(statistical_measure_t) :: rho_average

    call handle_cli()

    call init_main()

    !$omp parallel do schedule(dynamic) private(dynamics_state) firstprivate(time_control)
    do i_sample = 1, n_samples
        call init_dynamics(i_sample + (n_samples + 1))

        deallocate(dynamics_state)
    end do
    !$omp end parallel do

contains

    subroutine handle_cli()
        ! Default parameters
        rnd_seed = 9342342
        n_samples = 10
        edges_file = 'example.edgelist'
        par_b = 0.5_dp
        par_theta = 0.5_dp
        initial_infected_fraction = 1.0_dp
        beta_1 = 0.5_dp
        tmax = 100.0_dp
        time_scale = 'uniform'
        algorithm = 'HB_OGA'
        sampler_choice = 'rejection_maxheap'

        ! Here you can add code to handle command line arguments to override defaults
    end subroutine

    subroutine init_main()

        call time_control%init(time_scale)
        call time_average%init(time_control%get_max_array_size(real(tmax,dp)))
        call rho_average%init(time_control%get_max_array_size(real(tmax,dp)))

        ! Read the edgelist
        call read_network(net, edges_file)

        ! Clear and check the network (always necessary)
        call net%clear_and_check_all(min_order=1)

        ! Set the dynamical parameters
        call set_dyn_params(net, dyn_params, par_b, par_theta)

    end subroutine init_main

    subroutine init_dynamics(extra_seed)
        integer(kind=i4), intent(in), optional :: extra_seed
        integer(kind=i4) :: seed_offset

        if (present(extra_seed)) then
            seed_offset = extra_seed
        else
            seed_offset = 0
        end if

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

    subroutine set_dyn_params(net, dyn_params, par_b, par_theta)
        class(network_t), intent(in) :: net
        type(dyn_parameters_t), intent(inout) :: dyn_params
        real(kind=dp), intent(in) :: par_b, par_theta

        integer(kind=i4) :: edge_order

        ! Set parameters
        call dyn_params%init(net)

        ! Fill the beta and theta parameters
        do edge_order = 1, net%max_order
            dyn_params%beta(edge_order) = 1.0_dp + par_b * (edge_order -1)
            dyn_params%theta(edge_order) = ceiling(1.0_dp + (edge_order-1) * par_theta)
        end do

    end subroutine set_dyn_params

    subroutine read_network(net, edges_filename)
        type(network_t), intent(inout) :: net
        character(len=*), intent(in) :: edges_filename

        call network_import(net, trim(adjustl(edges_filename)))

    end subroutine read_network

end program simple_network_example_p
