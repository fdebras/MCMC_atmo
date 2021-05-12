import priors_atm


def make_prior(args):
    """Make a FullPrior for Jupiter for given command line options `args`"""
    prior_list = [
        priors_atm.SimplexUniformT_Priors(T_eq_min=800.0,T_eq_max=1500.0),#kappa_min=-5.0, kappa_max=-2.0, gamma_min=-2.0,gamma_max=2.0,T_int_min=100.0,T_int_max=300.0, \
    # T_eq_min=800.0,T_eq_max=1500.0),
        priors_atm.SimplexUniformAbund_Priors(H2O_min=-8.0, H2O_max=0.0,CO_min=-6.0, CO_max=0.0),
    # , CO2_min=-8.0, CO2_max=-2.0,CO_min=-8.0, CO_max=-2.0),
        priors_atm.SimplexUniformTransit_Priors(Kp_min=0.0,Kp_max=300.0, \
    Vsys_min=-20.0,Vsys_max=20.0)
        
    ]
#     if args.enable_diff_rotation:
#         prior_list.append(priors.DifferentialRotationPrior(omega_surface_absmax=0.2))

    return priors_atm.FullPrior(prior_list)
