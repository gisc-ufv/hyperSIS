module hyperSIS_program_common_mod
    use hyperSIS_kinds_mod, only: dp, i4
    use hyperSIS_network_mod, only: network_t
    use hyperSIS_network_io_mod, only: network_import
    use hyperSIS_dynamics_mod, only: dyn_parameters_t, net_state_base_t

    use rndgen_mod
    implicit none
    private

    interface
        subroutine proc_net_state_gen(net, state, gen)
            import :: net_state_base_t, rndgen, network_t
            class(network_t), intent(in) :: net
            class(net_state_base_t), intent(inout) :: state
            class(rndgen), intent(inout) :: gen
        end subroutine proc_net_state_gen
    end interface

    public :: proc_net_state_gen
    public :: check_qs_method
    public :: read_network
    public :: set_dyn_params
    public :: set_initial_number_of_infected_nodes
    public :: get_path_prefix

contains

    !> QS reactivation (randomly chosen node)
    subroutine reactivate_random_node(net, state, gen)
        class(network_t), intent(in) :: net
        class(net_state_base_t), intent(inout) :: state
        class(rndgen), intent(inout) :: gen

        integer(kind=i4) :: node_id

        if (state%get_num_infected() == 0) then
            node_id = gen%int(1, net%num_nodes)
            call state%add_infected(net, node_id)
        end if
    end subroutine reactivate_random_node

    !> Do nothing after
    subroutine do_nothing(net, state, gen)
        class(network_t), intent(in) :: net
        class(net_state_base_t), intent(inout) :: state
        class(rndgen), intent(inout) :: gen
    end subroutine do_nothing

    !> Check which method to use after dynamics step
    subroutine check_qs_method(after_dynamics_step, input_use_qs)
        procedure(proc_net_state_gen), pointer :: after_dynamics_step
        logical, intent(in) :: input_use_qs

        if (input_use_qs) then
            after_dynamics_step => reactivate_random_node
        else
            after_dynamics_step => do_nothing
        end if
    end subroutine

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

    subroutine set_initial_number_of_infected_nodes(net, inf_fraction, initial_number)
        class(network_t), intent(in) :: net
        real(kind=dp), intent(inout) :: inf_fraction
        integer(kind=i4), intent(inout) :: initial_number

        ! Check initial number of infected nodes
        if (initial_number > 0) then
            inf_fraction = real(initial_number,dp) / real(net%num_nodes,dp)
            if (inf_fraction < 0.0_dp .or. inf_fraction > 1.0_dp) error stop 'Initial number of infected nodes must be between 0 and network size'
        end if

    end subroutine set_initial_number_of_infected_nodes

    !> Build the prefix and remove old files
    function get_path_prefix(tmp_prefix, remove_files) result(prefix)
        character(len=*), intent(in) :: tmp_prefix
        logical, intent(in) :: remove_files
        character(len=:), allocatable :: prefix

        prefix = trim(adjustl(tmp_prefix))

        if (.not. prefix(len(prefix):len(prefix)) == '/') then
            ! Append a trailing underscore if it is not a directory
            prefix = prefix//'_'
        end if

        ! Remove old files with the same prefix
        call system('mkdir -pv '//prefix//'basedir/')
        if (remove_files) call system('rm -vf '//prefix//'*time*.dat')
    end function get_path_prefix

end module hyperSIS_program_common_mod
