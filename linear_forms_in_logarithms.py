def linear_forms_in_logarithms_bound(k, A, B, m=30, C1=22.8):
    """
    INPUT:
        k: a positive integer with k=2mod(4)
        A: a rational number
        B: a rational number
        m: the constant in Laurent's paper in Table 1
        C1: the constant in Laurent's paper in Table 1

    OUTPUT:
        An upper bound of n using linear forms in logarithms.
    """

    if k % 4 != 2:
        raise ValueError(f"k is not equal to 2 mod(4)!")

    a2 = A/2**((k/2) - 2)
    logA2 = max(log(a2.height()), 1)

    bound = 5

    ##### Case logA1 = 1 #####

    ### Max is m ###




    return logA2

