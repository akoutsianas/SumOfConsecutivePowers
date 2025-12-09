from sage.parallel.use_fork import p_iter_fork

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


class SumOfConsecutivePowersModularMethod:

    def __init__(self, k):
        self.k = k
        self.info = None
        self._x = polygen(QQ, 'x')
        self.fk = self._compute_fk()

    def _compute_fk(self):
        x = polygen(QQ, 'x')
        g = ((self._x - 1)**self.k + (self._x + 1)**self.k)/2
        f = self._x.parent(g/(self._x**2 + 1))
        f = sum([fi*x**i for i, fi in enumerate(f.coefficients())])
        return f

    def eliminate_newforms_method_1(self, primes_bound = 50):
        """
        INPUT
            - primes_bound: an upper bound of the primes we use in the elimination step

        OUTPUT:
            Saves in info variable a dictionary with the rational newforms over this elimination step fails and a list
            of small primes the modular method can not eliminate for the remaining newforms.
        """

        if self.k % 4 != 2:
            raise ValueError(f"k is not equal to 2 mod(4)!")

        Ex = lambda x: EllipticCurve([0, 2*x, 0, x**2 + 1, 0])
        values_d_d1 = VALUES_d_d1[self.k]
        info = {}
        for d in values_d_d1.keys():
            info[d] = {'failed_newforms': [], 'Bd': [2, 3, 5]}
            D = prod(ZZ(d).prime_factors())
            newfs = Newforms(2**7 * D, names="theta")
            print(f"We apply elimination step 1 for d={d}.")
            for newf in newfs:
                is_rational = 1 if newf.base_ring().degree == 1 else 0
                print(f"Elimination step for newform f={newf}")
                Bqs = []
                for q in prime_range(primes_bound):
                    if (2*D) % q != 0:
                        Bq = q
                        aqf = newf[q]
                        for x0 in range(q):
                            if (x0**2 + 1) % q == 0:
                                Bq *= aqf**2 - (q + 1)**2
                            else:
                                aEx0 = q + 1 - Ex(x0).reduction(q).order()
                                Bq *= aqf - aEx0
                        if is_rational:
                            Bqs.append(Bq)
                        else:
                            Bqs.append(Bq.norm())
                if gcd(Bqs) != 0:
                    for q in ZZ(gcd(Bqs)).prime_factors():
                        if q not in info[d]['Bd']:
                            info[d]['Bd'].append(q)
                else:
                    info[d]['failed_newforms'].append(newf)
        self.info = info

    def eliminate_newforms_method_2(self, bound_n, bound_t=50, lower_bound_n=7, ncpus=1):
        """
        INPUT:
            - bound_n: an upper bound of n
            - t_bound: a bound of t such that l = t*n + 1

        OUTPUT:
            A list of prime between 7 and bound_n that we are not able to eliminate using Proposition 5.3
        """
        if self.info is None:
            raise ValueError(f"You have to apply the 1st elimination method!")

        fork_iterator = p_iter_fork(ncpus=ncpus)
        fk = self.fk(self._x**2)
        fk /= 2**(self.k-2)
        problematic_n = []
        inputs = [([n], {'bound_t': bound_t}) for n in prime_range(lower_bound_n, bound_n + 1)]
        results = list(fork_iterator(self._eliminate_n, inputs))
        for result in  results:
            success_n = result[1]
            n = result[0][0][0]
            if not success_n:
                problematic_n.append(result[0][0][0])
            print(f"eliminate_n: {success_n}-{n}")

        return problematic_n

    def _eliminate_n(self, n, bound_t=50):
        print(f"Exponent n: {n}")
        for d in VALUES_d_d1[self.k].keys():
            for d1 in VALUES_d_d1[self.k][d]:
                if not self._eliminate_n_newforms(n, d, d1, bound_t=bound_t):
                    return False
        return True

    def _eliminate_n_newforms(self, n, d, d1, bound_t=50):
        Ex = lambda x: EllipticCurve([0, 2 * x, 0, x ** 2 + 1, 0])
        not_eliminated_newforms = self.info[d]['failed_newforms'].copy()
        for t in range(2, bound_t + 1):
            l = ZZ(t) * n + 1
            if (l in Primes()) and (ZZ(self.k / 2) % l != 0):
                for newf in self.info[d]['failed_newforms']:
                    if newf not in not_eliminated_newforms:
                        continue
                    suitable_newf_l = True
                    Fl = FiniteField(l)
                    fkbar = self.fk.change_ring(Fl)
                    y = polygen(Fl, 'y')
                    t_unit_roots = [r[0] for r in (y ** t - 1).roots()]
                    d2 = 1 / (Fl(d) ** 2 * Fl(d1))
                    alf = newf[l]
                    for zt in t_unit_roots:
                        x0s = [r[0] for r in (y ** 2 + 1 - 2 * Fl(d) * Fl(d1) * zt).roots()]
                        for x0 in x0s:
                            y2 = Fl(fkbar(x0) / (Fl(d) * Fl(d2)))
                            if y2 in t_unit_roots or y2.is_zero():
                                aEx0 = l + 1 - Ex(x0).order()
                                diff = (aEx0 - alf).norm()
                                if diff % n == 0:
                                    suitable_newf_l = False
                                    break
                        if not suitable_newf_l:
                            break
                    if l % 4 == 1:
                        diff = (4 - alf ** 2).norm()
                        if diff % n == 0:
                            suitable_newf_l = False
                    if suitable_newf_l:
                        not_eliminated_newforms.remove(newf)
                if len(not_eliminated_newforms) == 0:
                    return True
        return False
