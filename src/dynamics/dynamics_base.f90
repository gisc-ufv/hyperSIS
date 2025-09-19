module hyperSIS_dynamics_base_mod
    use hyperSIS_kinds_mod
    use hyperSIS_network_mod
    use datastructs_mod, only: dynamical_list_t
    implicit none
    private

    type, abstract :: state_compartment_base_t
        integer(kind=i4) :: num_nodes
        type(dynamical_list_t) :: nodes ! nodes in that state
        logical, allocatable :: is_edge_active(:) ! if the edge contains nodes in that state and is active
    end type

    type :: dyn_parameters_t
        real(kind=dp) :: alpha ! Infection rate for each node
        real(kind=dp), allocatable, dimension(:) :: beta ! Infection rate for each edge of order m
        integer(kind=i4), allocatable, dimension(:) :: theta ! Threshold for active edges of order m
        real(kind=dp) :: beta_scale
    contains
        procedure :: init => dyn_parameters_init
        procedure :: max_num_susceptible => dyn_parameters_max_num_susceptible
    end type

    type, abstract :: net_state_base_t
        integer(kind=i2), allocatable, dimension(:) :: node_state ! State of each node
        type(dyn_parameters_t) :: params
        real(kind=dp) :: time ! Time of the state
        real(kind=dp) :: total_rate, dt
    contains
        ! base procedures (to be overridden by extended types)
        procedure(net_state_base_init), deferred :: init
        procedure(net_state_base_add_infected), deferred :: add_infected
        procedure(net_state_base_remove_infected), deferred :: remove_infected
        procedure(net_state_base_dynamics_init), deferred :: dynamics_init
        procedure(net_state_base_dynamics_update_dt), deferred :: dynamics_update_dt
        procedure :: just_update_dt => dynamics_just_update_dt
        procedure(net_state_base_dynamics_step), deferred :: dynamics_step

        ! procedure to initialize the infected nodes (no need to override)
        procedure, private, non_overridable :: net_state_init_config_node
        procedure, private, non_overridable :: net_state_init_config_list_of_nodes
        procedure, private, non_overridable :: net_state_init_random_fraction_of_nodes
        generic :: init_config => net_state_init_config_node, net_state_init_config_list_of_nodes, net_state_init_random_fraction_of_nodes

        ! procedure to export the nodes states
        procedure :: export_nodes_states => net_state_export_nodes_states

        procedure(get_num_infected_interface), deferred :: get_num_infected
    end type

    abstract interface
        function get_num_infected_interface(this) result(n)
            import :: net_state_base_t, i4
            class(net_state_base_t), intent(in) :: this
            integer(i4) :: n
        end function
        ! base procedures
        subroutine net_state_base_init(this, net, params, sampler_choice)
            import :: net_state_base_t, network_t, dyn_parameters_t
            ! It will initialize the net_state with all nodes in empty state
            class(net_state_base_t) :: this
            type(network_t), intent(in) :: net
            class(dyn_parameters_t), intent(in) :: params
            character(len=*), intent(in) :: sampler_choice
        end subroutine

        subroutine net_state_base_dynamics_init(this, net)
            import :: net_state_base_t, network_t
            class(net_state_base_t) :: this
            class(network_t), intent(in) :: net
        end subroutine

        subroutine net_state_base_dynamics_step(this, net, gen)
            use rndgen_mod
            import :: net_state_base_t, network_t
            class(net_state_base_t) :: this
            class(network_t), intent(in) :: net
            type(rndgen) :: gen
        end subroutine

        function net_state_base_dynamics_update_dt(this, net, gen) result(res)
            ! will return .true. if the time step was updated
            ! will return .false. if the time step was not updated (for example, if there are no infected nodes)
            use rndgen_mod
            import :: net_state_base_t, network_t
            class(net_state_base_t) :: this
            class(network_t), intent(in) :: net
            type(rndgen) :: gen
            logical :: res
        end function

        subroutine net_state_base_add_infected(this, net, node_id)
            import :: net_state_base_t, network_t, i4
            class(net_state_base_t) :: this
            type(network_t), intent(in) :: net
            integer(kind=i4), intent(in) :: node_id
        end subroutine

        subroutine net_state_base_remove_infected(this, net, node_id, node_pos, new_state)
            import :: net_state_base_t, network_t, i4, i2
            class(net_state_base_t) :: this
            type(network_t), intent(in) :: net
            integer(kind=i4), intent(in) :: node_id, node_pos
            integer(kind=i2), intent(in) :: new_state
        end subroutine
        !/ base procedures
    end interface

    public :: net_state_base_t, dyn_parameters_t, state_compartment_base_t

contains

    ! just update dt (without returning a logical)
    subroutine dynamics_just_update_dt(this, net, gen)
        use rndgen_mod
        class(net_state_base_t) :: this
        class(network_t), intent(in) :: net
        type(rndgen) :: gen
        logical :: res

        res = this%dynamics_update_dt(net, gen)
    end subroutine

    ! subroutines to initialize the infected nodes
    subroutine net_state_init_config_node(this, net, params, sampler_choice, node_id)
        ! Initialize the state with a random node
        class(net_state_base_t) :: this
        type(network_t), intent(in) :: net
        integer(kind=i4), intent(in) :: node_id
        class(dyn_parameters_t), intent(in) :: params
        character(len=*), intent(in) :: sampler_choice

        call net_state_init_config_list_of_nodes(this, net, params, sampler_choice, [node_id])

    end subroutine

    subroutine net_state_init_config_list_of_nodes(this, net, params, sampler_choice, nodes)
        ! Initialize the state with a random node
        class(net_state_base_t) :: this
        type(network_t), intent(in) :: net
        integer(kind=i4), intent(in) :: nodes(:)
        class(dyn_parameters_t), intent(in) :: params
        integer(kind=i4) :: node_pos
        character(len=*), intent(in) :: sampler_choice

        call this%init(net, params, sampler_choice)
        do node_pos = 1, size(nodes)
            call this%add_infected(net, nodes(node_pos))
        end do

    end subroutine

    subroutine net_state_init_random_fraction_of_nodes(this, net, gen, params, sampler_choice, fraction)
        use rndgen_mod
        ! Initialize the state with a random fraction of nodes
        class(net_state_base_t) :: this
        type(network_t), intent(in) :: net
        type(rndgen), intent(in) :: gen
        class(dyn_parameters_t), intent(in) :: params
        real(kind=dp), intent(in) :: fraction
        integer(kind=i4) :: max_num_nodes, node_id, cnt_num_nodes
        character(len=*), intent(in) :: sampler_choice

        call this%init(net, params, sampler_choice)

        max_num_nodes = int(fraction * net%num_nodes)
        cnt_num_nodes = 0

        if (max_num_nodes == net%num_nodes) then
            do node_id = 1, net%num_nodes
                call add_infected(node_id)
            end do
        end if

        do while (cnt_num_nodes < max_num_nodes)
            ! Select a random node
            call add_infected(gen%int(1, net%num_nodes))
        end do

        if (this%get_num_infected() == 0) error stop 'No infected nodes found'

    contains

        subroutine add_infected(node_id)
            integer(kind=i4), intent(in) :: node_id

            if (this%node_state(node_id) == 0_i2) then
                ! If the node is not infected, infect it
                call this%add_infected(net, node_id)
                cnt_num_nodes = cnt_num_nodes + 1
            end if
        end subroutine

    end subroutine
    !/ subroutines to initialize the infected nodes

    ! subroutines to export the states
    subroutine net_state_export_nodes_states(this, net, filename)
        ! Export the nodes states to a file
        class(net_state_base_t) :: this
        type(network_t), intent(in) :: net
        character(len=*), intent(in) :: filename
        integer(kind=i4) :: funit, i

        open(newunit=funit, file=filename, status='replace', action='write')
        do i = 1, net%num_nodes
            write(funit, fmt_general) i, this%node_state(i)
        end do
        close(funit)

    end subroutine

    ! subroutines for the dynamical parameters
    subroutine dyn_parameters_init(this, net)
        class(dyn_parameters_t) :: this
        type(network_t) :: net

        allocate(this%beta(net%max_order))
        allocate(this%theta(net%max_order))

        ! Set default scale
        this%beta_scale = 1.0_dp

        ! Set default values
        this%alpha = 1.0_dp
        this%beta = 1.0_dp
        this%theta = 1

    end subroutine

    function dyn_parameters_max_num_susceptible(this, edge_order) result(res)
        ! Function to return the maximum number of susceptible nodes in an edge of order m
        class(dyn_parameters_t) :: this
        integer(kind=i4), intent(in) :: edge_order
        integer(kind=i4) :: res

        res = (edge_order + 1 - this%theta(edge_order))
    end function
    !/ subroutines for the dynamical parameters

end module
