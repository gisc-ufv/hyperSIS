module hyperSIS_dynamics_mod
    use hyperSIS_dynamics_chooser_mod, only: net_state_choose, dynamics_choices => alg_choices
    use hyperSIS_dynamics_base_mod, only: dyn_parameters_t, net_state_base_t
    implicit none
    private

    public :: net_state_base_t, dyn_parameters_t, net_state_choose, dynamics_choices
end module
