module hyperSIS_dynamics_HB_OGA_mod
    use hyperSIS_kinds_mod
    use hyperSIS_network_mod
    use hyperSIS_dynamics_base_mod
    use datastructs_mod, only: dynamical_list_t
    use datastructs_mod, only: sampler_base_t
    use datastructs_mod, only: choose_sampler
    use datastructs_mod, only: log_unit, log_write, LOGGER_OK, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG
    implicit none
    private

    character(len=*), parameter :: ALGORITHM_NAME = 'HB_OGA'

    type, extends(state_compartment_base_t) :: state_compartment_t
        integer(kind=i4), allocatable :: num_infected_nodes_per_edge(:) ! number of nodes in that state in that edge
    end type

    type, extends(net_state_base_t) :: net_state_t
        type(state_compartment_t) :: infected ! Infected nodes, base type
        ! HB-OGA specific variables
        real(kind=dp) :: total_infection_attempt_rate ! Rate of infection attempts
        real(kind=dp) :: total_healing_rate ! Rate of healing
        class(sampler_base_t), allocatable :: infection_attempts_sampler ! Sampler for the infection attempts
    contains
        procedure :: init => net_state_init
        procedure :: add_infected => net_state_add_infected
        procedure :: remove_infected => net_state_remove_infected
        ! HB-OGA will be the default algorithm
        procedure :: dynamics_init => net_state_dynamics_init
        procedure :: dynamics_update_dt => net_state_dynamics_update_dt
        procedure :: dynamics_step => net_state_dynamics_step
        procedure :: calculate_rates => net_state_calculate_rates
        procedure :: print_debug_quantities => net_state_print_debug_quantities
        procedure :: remove_possibly_active_edge => net_state_remove_possibly_active_edge_from_list
        procedure :: activate_edge => net_state_activate_edge

        procedure :: export_edges_states => net_state_export_edges_states

        procedure :: get_num_infected => net_state_get_num_infected
    end type

    public :: net_state_t

contains

    function net_state_get_num_infected(this) result(n)
        ! Get the number of infected nodes
        class(net_state_t), intent(in) :: this
        integer(kind=i4) :: n

        n = this%infected%num_nodes

    end function net_state_get_num_infected

    subroutine net_state_export_edges_states(this, net, filename)
        ! Export the edges states to a file
        class(net_state_t) :: this
        type(network_t), intent(in) :: net
        character(len=*), intent(in) :: filename
        integer(kind=i4) :: funit, i

        open(newunit=funit, file=filename, status='replace', action='write')
        write(funit, fmt_general) "# time = ", this%time
        do i = 1, net%num_edges
            if (this%infected%is_edge_active(i)) write(funit, fmt_general) i, merge(1, 0, this%infected%is_edge_active(i))
        end do
        close(funit)

    end subroutine

    subroutine net_state_init(this, net, params, sampler_choice)
        ! It will initialize the net_state with all nodes in empty state
        class(net_state_t) :: this
        type(network_t), intent(in) :: net
        class(dyn_parameters_t), intent(in) :: params
        character(len=*), intent(in) :: sampler_choice

        call log_write(LOG_INFO, 'Dynamics: '//trim(adjustl(ALGORITHM_NAME)), .false.)
        call log_write(LOG_INFO, 'Sampler: '//trim(adjustl(sampler_choice)))

        ! allocate sampler
        call choose_sampler(this%infection_attempts_sampler, sampler_choice)

        this%params = params
        this%infected%num_nodes = 0

        allocate(this%node_state(net%num_nodes))
        allocate(this%infected%num_infected_nodes_per_edge(net%num_edges))
        allocate(this%infected%is_edge_active(net%num_edges))

        this%node_state = 0
        this%infected%num_infected_nodes_per_edge = 0
        this%infected%is_edge_active = .false.

        ! initialize dynamical lists
        call this%infected%nodes%init(net%num_nodes)

        ! initialize the infection attempts sampler
        select case (trim(adjustl(sampler_choice)))
          case ('rejection_maxheap_composition')
            block
                real(kind=dp) :: min_weight, max_weight, weight
                integer(kind=i4) :: edge_order

                min_weight = huge(min_weight)
                max_weight = 0.0_dp

                do edge_order = 1, net%max_order
                    weight = this%params%beta(edge_order) * this%params%max_num_susceptible(edge_order)
                    if (weight < min_weight) min_weight = weight
                    if (weight > max_weight) max_weight = weight
                end do

                call this%infection_attempts_sampler%init(net%num_edges, min_weight, max_weight)
            end block
          case ('rejection_two_classes', 'rejection_maxheap_two_classes')
            block
                real(kind=dp) :: threshold, threshold_sum, threshold_std, weight
                integer(kind=i4) :: edge_order

                threshold_sum = 0.0_dp
                threshold_std = 0.0_dp
                do edge_order = 1, net%max_order
                    weight = this%params%beta(edge_order) * this%params%max_num_susceptible(edge_order)
                    threshold_sum = threshold_sum + weight
                    threshold_std = threshold_std + weight**2
                end do

                threshold_std = sqrt(threshold_std / net%max_order - (threshold_sum / net%max_order)**2)

                threshold = threshold_sum / net%max_order + threshold_std

                write(*, fmt_general) 'Initializing with threshold: ', threshold

                call this%infection_attempts_sampler%init(net%num_edges, threshold)
            end block
          case default
            call this%infection_attempts_sampler%init(net%num_edges)
        end select

    end subroutine

    subroutine net_state_dynamics_init(this, net)
        ! Initialize the HB-OGA dynamics
        class(net_state_t) :: this
        class(network_t), intent(in) :: net ! compatibility with the other algorithms

        ! Set the time to zero
        this%time = 0.0_dp

        call this%calculate_rates(net)

    end subroutine

    subroutine net_state_calculate_rates(this, net)
        ! Calculate the rates of infection and healing
        class(net_state_t) :: this
        class(network_t), intent(in) :: net

        ! Initial healing rate is the number of infected nodes times the healing rate
        this%total_healing_rate = this%infected%nodes%n_used * this%params%alpha

        ! Loop over the active edges and calculate the total infection attempt rate
        this%total_infection_attempt_rate = this%params%beta_scale * this%infection_attempts_sampler%sum()

    end subroutine net_state_calculate_rates

    function net_state_dynamics_update_dt(this, net, gen) result(res)
        use rndgen_mod
        class(net_state_t) :: this
        class(network_t), intent(in) :: net
        type(rndgen) :: gen
        logical :: res

        res = .true.

        if (this%infected%num_nodes == 0) then
            ! no infected nodes, nothing to do
            print*, 'No infected nodes, nothing to do.'
            res = .false.
            return
        end if

        call this%calculate_rates(net)

        this%total_rate = this%total_infection_attempt_rate + this%total_healing_rate

        ! Calculate the time step
        this%dt = -log(1.0_dp - gen%rnd()) / this%total_rate

        ! Update the time
        this%time = this%time + this%dt

    end function

    subroutine net_state_dynamics_step(this, net, gen)
        use rndgen_mod
        ! Perform a step of the HB-OGA dynamics
        class(net_state_t) :: this
        class(network_t), intent(in) :: net
        type(rndgen) :: gen
        integer(kind=i4) :: node_pos, node_id
        integer(kind=i4) :: edge_id, edge_order, edge_num_susceptible

        if (gen%rnd() < this%total_healing_rate/this%total_rate) then
            ! healing will happen
            ! we select a random infected node

            node_pos = gen%int(1, this%infected%nodes%n_used)
            node_id = this%infected%nodes%list(node_pos)

            ! remove it from the infected state and move to susceptible state
            call this%remove_infected(net, node_id, node_pos, 0_i2)
        else
            ! infection will happen
            ! for that, we select a random active edge using rejection sampling

            edge_id = this%infection_attempts_sampler%sample(gen)
            edge_order = net%edges(edge_id)%order

            ! if not active, we only remove it from the list
            if (.not. this%infected%is_edge_active(edge_id)) then
                ! we only remove that edge from the list (phantom process)
                call this%remove_possibly_active_edge(edge_id, edge_order)
                return ! nothing to do, we just remove the edge from the list
            end if

            ! here, the edge is indeed active, so we will try to infect a node in that edge

            !print*, 'accepted edge', edge_id, 'of order', edge_order, 'with', this%infected%num_infected_nodes_per_edge(edge_id), 'nodes in that state'

            edge_num_susceptible = edge_order + 1 - this%infected%num_infected_nodes_per_edge(edge_id)

            !print*, edge_num_susceptible, this%params%max_num_susceptible(edge_order)

            ! now we know the active edge, and see if we will indeed infect a node
            if (gen%rnd() < 1.0_dp*(edge_num_susceptible) / this%params%max_num_susceptible(edge_order) ) then
                ! we will infect a node in that edge, using a list of susceptible individuals
                block
                    integer(kind=i4) :: list_of_sus(edge_num_susceptible)
                    integer(kind=i4) :: sus_pos

                    sus_pos = 0

                    do node_pos = 1, net%edges(edge_id)%order + 1
                        node_id = net%edges(edge_id)%nodes(node_pos)
                        if (this%node_state(node_id) == 0_i2) then
                            sus_pos = sus_pos + 1
                            list_of_sus(sus_pos) = node_id
                        end if
                    end do

                    ! DEBUG
                    !if (list_of_sus%n_used /= edge_num_susceptible) stop 'eita'

                    node_pos = gen%int(1, edge_num_susceptible)
                    node_id = list_of_sus(node_pos)

                    call this%add_infected(net, node_id)
                end block
            else
                ! otherwise, we do nothing! Phantom process :)
            end if

        end if

    end subroutine

    ! @@@@@@@@@@@@@@@@ HB_OGA - auxiliary functions @@@@@@@@@@@@@@@@@@
    subroutine check_infected_edge_condition(this, net, edge_id, edge_order)
        ! Subroutine to check if the edge is active in the infected compartment
        class(net_state_t) :: this
        type(network_t), intent(in) :: net
        integer(kind=i4), intent(in) :: edge_id, edge_order

        if ( &
            (this%infected%num_infected_nodes_per_edge(edge_id) >= this%params%theta(edge_order)) & ! it has to have at least theta nodes in that state
            .and. (this%infected%num_infected_nodes_per_edge(edge_id) /= net%edges(edge_id)%order + 1) & ! it has to have at least one node in another state
            ) then
            ! it will be active
            if (.not. this%infected%is_edge_active(edge_id)) then
                !print(fmt_general), 'Activating edge', edge_id, 'of order', edge_order, 'with', this%infected%num_infected_nodes_per_edge(edge_id), 'nodes in that state'
                call this%activate_edge(edge_id, edge_order)
            endif
            ! otherwise, it is already active
        else
            ! only change the edge status, do nothing in the list
            this%infected%is_edge_active(edge_id) = .false.
        end if

    end subroutine

    subroutine net_state_add_infected(this, net, node_id)
        ! Subroutine to add infected node in the list
        ! We assume that it was not infected before, be careful!
        class(net_state_t) :: this
        type(network_t), intent(in) :: net
        integer(kind=i4), intent(in) :: node_id
        integer(kind=i4) :: edge_pos, edge_id, edge_order

        ! DEBUG
        !if (this%node_state(node_id) == 1_i2) stop 'Node already infected!'

        this%node_state(node_id) = 1_i2 ! add to that state

        ! add node to the list
        this%infected%num_nodes = this%infected%num_nodes + 1
        call this%infected%nodes%add(node_id)

        ! fill the edges list
        ! if the edge is not active, add it. Otherwise, keep as is
        do edge_pos = 1, net%nodes(node_id)%degree
            edge_id = net%nodes(node_id)%edges(edge_pos)
            edge_order = net%edges(edge_id)%order
            ! add one node to that state in that edge
            this%infected%num_infected_nodes_per_edge(edge_id) = this%infected%num_infected_nodes_per_edge(edge_id) + 1
            ! if the edge is not active, it can become active
            ! else, if the edge is active, it can become inactive
            call check_infected_edge_condition(this, net, edge_id, edge_order)
        end do

    end subroutine

    subroutine net_state_remove_infected(this, net, node_id, node_pos, new_state)
        ! Subroutine to remove infected node in the list
        ! CAUTION: WILL NOT UPDATE THE ACTIVE EDGES LIST
        class(net_state_t) :: this
        type(network_t), intent(in) :: net
        integer(kind=i4), intent(in) :: node_id, node_pos
        integer(kind=i4) :: edge_pos, edge_id, edge_order
        integer(kind=i2), intent(in) :: new_state

        ! DEBUG
        !if (this%node_state(node_id) /= 1_i2) stop 'Node was not infected!'
        !if (this%infected%nodes%list(node_pos) /= node_id) stop 'Node is not the correct one!'

        this%node_state(node_id) = new_state ! add to that state

        ! remove node from the list
        this%infected%num_nodes = this%infected%num_nodes - 1
        call this%infected%nodes%remove(node_pos)

        ! update the status of the edges
        do edge_pos = 1, net%nodes(node_id)%degree
            edge_id = net%nodes(node_id)%edges(edge_pos)
            edge_order = net%edges(edge_id)%order
            this%infected%num_infected_nodes_per_edge(edge_id) = this%infected%num_infected_nodes_per_edge(edge_id) - 1
            ! the number of active edges can be reduced or increased
            ! we check the condition
            call check_infected_edge_condition(this, net, edge_id, edge_order)
        end do

    end subroutine

    subroutine net_state_print_debug_quantities(this)
        ! Print debug quantities of the network state
        class(net_state_t) :: this

        print(fmt_general), 'Number of infected nodes:', this%infected%num_nodes
        print(fmt_general), 'Number of active edges:', count(this%infected%is_edge_active)
        !print(fmt_general), 'Number of active edges in the list:', this%infected%possibly_active_edges%n_used
        !print(fmt_general), 'Maximum order of possibly active edges:', this%infected%possibly_active_max_order
        print(fmt_general), 'Total infection attempt rate:', this%total_infection_attempt_rate
        print(fmt_general), 'Total healing rate:', this%total_healing_rate
        print(fmt_general), 'Infection attempt probability:', &
            this%total_infection_attempt_rate / (this%total_infection_attempt_rate + this%total_healing_rate)
        print(fmt_general), ''

    end subroutine

    subroutine net_state_remove_possibly_active_edge_from_list(this, edge_id, edge_order)
        class(net_state_t) :: this
        integer(kind=i4), intent(in) :: edge_id, edge_order
        ! Remove the edge from the active edges list

        ! DEBUG
        !if (.not. this%infected%is_edge_in_weighted_sampler(edge_id)) stop 'Edge is not in the active edges list!'
        !if (this%infected%is_edge_active(edge_id)) stop 'Edge is active!'

        call this%infection_attempts_sampler%remove(edge_id)
    end subroutine

    subroutine net_state_activate_edge(this, edge_id, edge_order)
        ! Subroutine to activate the edge in the compartment
        class(net_state_t) :: this
        integer(kind=i4), intent(in) :: edge_id, edge_order

        ! DEBUG
        !if (this%infected%is_edge_active(edge_id)) stop 'Edge is already active!'

        this%infected%is_edge_active(edge_id) = .true.

        call this%infection_attempts_sampler%set_weight(edge_id, this%params%beta(edge_order) * this%params%max_num_susceptible(edge_order))
    end subroutine
    !/@@@@@@@@@@@@@@@@ HB_OGA @@@@@@@@@@@@@@@@@@

end module
