VALUES_d_d1 = {
    10: {5: [1, 1/ZZ(5)**2]},
    26: {13: [1, 1/ZZ(13)**2]},
    30: {5: [1, 1/ZZ(5)**2]},
    34: {17: [1, 1/ZZ(17)**2]},
    50: {5: [1, 1/ZZ(5)**2], 25: [1, 1/ZZ(5)**4]},
    58: {29: [1, 1/ZZ(29)**2]},
    70: {5: [1, 1/ZZ(5)**2]},
    74: {37: [1, 1/ZZ(37)**2]},
    78: {13: [1, 1/ZZ(13)**2]},
    82: {41: [1, 1/ZZ(41)**2]},
    90: {5: [1, 1/ZZ(5)**2]},
}


def bound_exponent_n(k, primes_bound = 50):
    """
    INPUT
        - k: an integer equivalent 2 mod 4
        - primes_bound: an upper bound of the primes we use in the elimination step

    OUTPUT:
        A dictionary with the rational newforms over this elimination step fails and a list of small primes the modular
        method can not eliminate for the remaining newforms.
    """

    if k % 4 != 2:
        raise ValueError(f"k is not equal to 2 mod(4)!")

    Ex = lambda x: EllipticCurve([0, 2*x, 0, x**2 + 1, 0])
    values_d_d1 = VALUES_d_d1[k]
    info = {}
    for d in values_d_d1.keys():
        info[d] = {'failed_newforms': [], 'Bd': [2, 3, 5]}
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
            # print(f"gcd(Bqs): {gcd(Bqs)}")
            if gcd(Bqs) != 0:
                # print(f"small_primes: {ZZ(gcd(Bqs)).prime_factors()}")
                for q in ZZ(gcd(Bqs)).prime_factors():
                    if q not in info[d]['Bd']:
                        info[d]['Bd'].append(q)
            else:
                info[d]['failed_newforms'].append(fnew)
            # break
    return info


def eliminate_for_give_n(k, info, bound_n, bound_t=50):
    """
    INPUT:
        - k: an integer equivalent 2 mod 4
        - info: the dictionary that we get from bound_exponent_n function
        - Bn: an upper bound of n
        - t_bound: a bound of t such that l = t*n + 1

    OUTPUT:
        A list of prime betwenn 7 and Bn that we are not able to eliminate all newforms
    """
    fk = compute_fk(k)
    fk = fk(fk.parent().gen()**2)
    print(f"fk: {fk}")
    fk /= 2**(k-2)
    Ex = lambda x: EllipticCurve([0, 2 * x, 0, x ** 2 + 1, 0])
    for n in prime_range(7, bound_n + 1):
        print(f"Eliminate n: {n}")
        eliminate_n = False
        for t in range(2, bound_t):
            l = ZZ(t)*n + 1
            if l in Primes() and ZZ(k/2) % l != 0:
                # print(f"l: {l}")
                suitable_l = True
                Fl = FiniteField(l)
                fkbar = fk.change_ring(Fl)
                y = polygen(Fl, 'y')
                t_unit_roots = [r[0] for r in (y**t - 1).roots()]
                for d in VALUES_d_d1[k].keys():
                    for d1 in VALUES_d_d1[k][d]:
                        d2 = 1/(Fl(d)**2 * Fl(d1))
                        for zt in t_unit_roots:
                            X0 = [r[0] for r in (y**2 + 1 - 2*Fl(d)*Fl(d1)*zt).roots()]
                            for x0 in X0:
                                if Fl(fkbar(x0) / (Fl(d) * Fl(d2))) in t_unit_roots:
                                    aEx0 = l + 1 - Ex(x0).order()
                                    for i, fnew in enumerate(info[d]['failed_newforms']):
                                        diff = aEx0 - fnew[l]
                                        if diff % n == 0:
                                            suitable_l = False
                                            break
                                        # print(f"diff-{i}: {diff}")
                # print(f"Suitable_l: {suitable_l}")
                if suitable_l:
                    eliminate_n = True
                    break
        print(f"eliminate_n: {eliminate_n}-{n}")



def compute_fk(k):
    x = polygen(QQ, 'x')
    g = ((x - 1)**k + (x + 1)**k)/2
    f = x.parent(g/(x**2 + 1))
    f = sum([fi*x**i for i, fi in enumerate(f.coefficients())])
    return f
