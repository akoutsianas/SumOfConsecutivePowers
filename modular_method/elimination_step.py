VALUES_d_d1 = {
    10: {
        5: [1, 1/5**2]
    },
    26: {
        13: [1, 1/13**2]
    }
}

def bound_exponent_n(k, primes_bound = 50):
    """
    INPUT
        - k: an integer equivalent 2 mod 4
        - primes_bound: an upper bound of the primes we use in the elimination step

    OUTPUT:
        - bol: True, if the modular method succeed, otherwise False
        - bound: an upper bound of the exponent n
    """

    if k % 4 != 2:
        raise ValueError(f"k is not 2 mod(4)!")

    Ex = lambda x: EllipticCurve([0, 2*x, 0, x**2 + 1, 0])
    k_odd = ZZ(k/2)
    values_d_d1 = VALUES_d_d1[k]
    for d in values_d_d1.keys():
        for d1 in values_d_d1[d]:
            D = prod(ZZ(d).prime_factors())
            newfs = Newforms(2**7 * D, names="theta")
            print(f"We apply elimination step for d={d}.")
            for fnew in newfs:
                is_rational = 1 if fnew.base_ring().degree == 1 else 0
                print(f"Elimination step for newform f={fnew}")
                Bqs = []
                for q in prime_range(primes_bound):
                    if (2*D) % q != 0:
                        # print(f"q: {q}")
                        Bq = q
                        aqf = fnew[q]
                        for x0 in range(q):
                            # print(f"x0: {x0}")
                            if (x0**2 + 1) % q == 0:
                                Bq *= aqf**2 - (q + 1)**2
                            else:
                                # print(f"Ex0: {Ex(x0)}")
                                aEx0 = q + 1 - Ex(x0).reduction(q).order()
                                # print(f"dif: {aqf - aEx0}")
                                Bq *= aqf - aEx0
                        # print(f"Bq: {Bq}")
                        if is_rational:
                            Bqs.append(Bq)
                        else:
                            Bqs.append(Bq.norm())
                        # print(f"Bqs: {Bqs}")
                        # print(f"------")
                # print(f"Bqs: {Bqs}")
                print(f"gcd(Bqs): {gcd(Bqs)}")
                if gcd(Bqs) != 0:
                    print(f"small_primes: {ZZ(gcd(Bqs)).prime_factors()}")
                # break

