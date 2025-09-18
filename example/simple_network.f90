program simple_network_example_p
    use hyperSIS_kinds_mod
    use hyperSIS_network_mod
    use hyperSIS_network_io_mod
    use hyperSIS_dynamics_mod

    use rndgen_mod
    implicit none

    type(network_t) :: net
    type(dyn_parameters_t) :: dyn_params
    class(net_state_base_t), allocatable :: dynamics_state

    type(rndgen) :: dyn_gen

    integer(kind=i4) :: i

    call dyn_gen%init(9342342) ! initialize the random generator with a seed

    ! Generate a simple network manually
    call generate_manual_network(net)

    ! OR read a network from an edgelist file
    ! (do not use both at the same time)
    !call read_network(net, 'example.edgelist')

    ! Clean and check the network (always necessary)
    call net%clear_and_check_all(min_order=1)

    ! Set the dynamical parameters
    call set_dyn_params(net, dyn_params, par_b=0.5_dp, par_theta=0.5_dp)

    ! select dynamics algorithm
    call net_state_choose(dynamics_state, 'HB_OGA')

    ! allocate and initializes the dynamics with a random configuration
    call dynamics_state%init_config(net=net, gen=dyn_gen, params=dyn_params, sampler_choice='rejection_maxheap', fraction=1.0_dp)

    ! Set the beta scale parameter
    ! by default, alpha = 1.0
    dynamics_state%params%beta_scale = 0.5_dp

    ! Initialize the dynamics
    call dynamics_state%dynamics_init(net)

    ! Perform 10000 steps of the dynamics
    do i = 1, 10000
        ! get Gillespie time
        if (.not. dynamics_state%dynamics_update_dt(net, dyn_gen)) exit ! exit if no infected nodes

        ! update the state
        call dynamics_state%dynamics_step(net, dyn_gen)

        ! Print the time and number of infected nodes
        write(*,fmt_general) 'Step ', i, ': time = ', dynamics_state%time, 1.0_dp * dynamics_state%get_num_infected() / net%num_nodes
    end do

contains

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

    subroutine generate_manual_network(net)
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

    end subroutine generate_manual_network

    subroutine read_network(net, edges_filename)
        type(network_t), intent(inout) :: net
        character(len=*), intent(in) :: edges_filename

        call network_import(net, trim(adjustl(edges_filename)))

    end subroutine read_network

end program simple_network_example_p