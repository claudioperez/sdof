#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
def parameters(
        radius=None, # spectral radius in the high-frequency limit, rho_inf
        family=None,
    ):

    # For the following to work, the arguments to ^
    # must both be bool
    #                       xor
    assert (alpha_m is None) ^ (radius is None)
    assert (alpha_f is None) ^ (radius is None)

    assert (alpha_m is None) !=  (radius is None)
    assert (alpha_f is None) !=  (radius is None)

    if gamma is None and beta is None:
        # 1) Given radius + family
        if radius is None:
            assert alpha_m is not None
            assert alpha_f is not None

        elif family == "HHT":
            alpha_f = (1-radius)/(1+radius) # (22)
            alpha_m = 0.0

        elif family == "WBZ":
            alpha_f = 0.0
            alpha_m = (1-radius)/(1+radius) # (23)

        elif family == "alpha":
            # Minimize low-freq dissipation by enforcing 
            # alpha_f = (alpha_m + 1)/3 :
            alpha_f = radius/(radius + 1)           # (25)
            alpha_m = (2*radius - 1)/(radius + 1)
        else:
            raise ValueError(f"Unknown family {family}")

        # 2) Given alpha_f + alpha_m

        gamma = 0.5 - alpha_m + alpha_f           # (17) for second-order accuracy
        beta  = 0.25*(1 - alpha_m + alpha_f)**2   # (20) 
    pass
