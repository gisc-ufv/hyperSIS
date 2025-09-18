module hyperSIS_dynamics_chooser_mod
    use hyperSIS_kinds_mod
    use hyperSIS_dynamics_base_mod, only: net_state_base_t
    implicit none
    private

    character(len=*), parameter :: alg_choices = 'HB_OGA,NB_OGA'

    public :: net_state_choose, alg_choices

contains

    subroutine net_state_choose(net_state, selected_algorithm)

        use hyperSIS_dynamics_HB_OGA_mod, only : HB_OGA_mod_net_state_t => net_state_t
        use hyperSIS_dynamics_NB_OGA_mod, only : NB_OGA_mod_net_state_t => net_state_t

        class(net_state_base_t), allocatable, intent(out) :: net_state
        character(len=*), intent(in) :: selected_algorithm

        select case (trim(adjustl(selected_algorithm)))
            case ('HB_OGA')
                allocate(HB_OGA_mod_net_state_t :: net_state)
            case ('NB_OGA')
                allocate(NB_OGA_mod_net_state_t :: net_state)
            case default
                error stop 'Unknown algorithm selected: '//trim(adjustl(selected_algorithm))
        end select
    end subroutine

end module
