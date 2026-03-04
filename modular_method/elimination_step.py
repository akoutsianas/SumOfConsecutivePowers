import time
from sage.all import polygen, ZZ, QQ, EllipticCurve, Newforms, prime_range, prod, gcd, Primes, FiniteField, Integer
from sage.parallel.use_fork import p_iter_fork
from sage.parallel.decorate import Parallel

from modular_method.bounds_info import pairs_info_n0_values


class SumOfConsecutivePowersModularMethod:

    def __init__(self, k):
        if k % 4 != 2:
            raise ValueError(f"k is not equal to 2 mod(4)!")
        self.k = k
        self.info = {}
        self._x = polygen(QQ, 'x')
        self.fk = self._compute_fk()
        self.pairs_d1_d2 = self._compute_pairs_d1_d2()

    def _compute_fk(self):
        g = ((self._x - 1)**self.k + (self._x + 1)**self.k)/2
        f = self._x.parent(g/(self._x**2 + 1))
        f = sum([fi * self._x**i for i, fi in enumerate(f.coefficients())])
        return f

    def _compute_pairs_d1_d2(self):
        divs = ZZ(self.k/2).divisors()
        divs = [d for d in divs if len([1 for p in d.prime_divisors() if p % 4 == 3]) == 0]
        pairs = []
        for d1 in divs:
            for d2 in divs:
                if gcd(d1, d2) == 1 and d1*d2 != 1:
                    pairs.append([d1, d2])
        return pairs

    def eliminate_newforms_method_1(self, primes_bound=50):
        for pair in self.pairs_d1_d2:
            d1 = pair[0]
            d2 = pair[1]
            self.info[f"pair_{d1}_{d2}"] = self._eliminate_newforms_method_1_d1_d2(d1, d2, primes_bound=primes_bound)

    def _eliminate_newforms_method_1_d1_d2(self, d1, d2, primes_bound=50):
        """
        INPUT
            - primes_bound: an upper bound of the primes we use in the elimination step

        OUTPUT:
            Saves in info variable a dictionary with the rational newforms over this elimination step fails and a list
            of small primes the modular method can not eliminate for the remaining newforms.
        """
        Ex = lambda x: EllipticCurve([0, 2*x, 0, x**2 + 1, 0])
        info = {'failed_newforms': [], 'Bd': [2, 3, 5]}
        D = prod(ZZ(d1*d2).prime_factors())
        newfs = Newforms(2**7 * D, names="theta")
        print(f"We apply elimination step 1 for d1={d1} and d2={d2}.")
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
                    if q not in info['Bd']:
                        info['Bd'].append(q)
            else:
                info['failed_newforms'].append(newf.abelian_variety().elliptic_curve())
        return info

    def eliminate_newforms_method_2(self, bound_t=150, lower_bound_n=7, ncpus=1):
        for pair_info in pairs_info_n0_values[self.k]:
            d1 = pair_info[0]
            d2 = pair_info[1]
            bound = pair_info[2]
            t0 = time.time()
            problematic_n = self.eliminate_newforms_method_2_d1_d2(d1, d2, bound, bound_t=bound_t,
                                                                   lower_bound_n=lower_bound_n, ncpus=ncpus)
            t1 = time.time()
            print(f"The problematic exponents for d1={d1} and d2=d2 are {problematic_n}. Time for computation={t1-t0} "
                  f"and exponents bound={bound}")

    def eliminate_newforms_method_2_d1_d2(self, d1, d2, exponents, bound_t=150, lower_bound_n=7, ncpus=1):
        """
        INPUT:
            - exponents: an upper bound of n or a list of exponents
            - bound_t: a bound of t such that l = t*n + 1
            - lower_bound_n: the lower bound of the range of primes
            - ncpus: the number of cpus we use in parallel

        OUTPUT:
            A list of prime between 7 and bound_n that we are not able to eliminate using Proposition 5.3
        """
        if self.info is None:
            raise ValueError(f"You have to apply the 1st elimination method!")
        print(f"Elimination method 2 for d1={d1} and d2={d2}.")

        fork_iterator = p_iter_fork(ncpus=ncpus)
        if isinstance(exponents, Integer):
            primes_range = prime_range(lower_bound_n, exponents + 1)
        elif isinstance(exponents, list):
            primes_range = exponents
        else:
            raise ValueError('Exponents is not an integer or a list!')
        inputs = [
            (
                [d1, d2, [n for j, n in enumerate(primes_range) if j % ncpus == i]],
                {'bound_t': bound_t}
            ) for i in range(ncpus)]
        results = list(fork_iterator(self._eliminate_list_of_n_d1_d2, inputs))
        problematic_n = []
        for result in results:
            for n in result[1]:
                problematic_n.append(n)
        return problematic_n

    def _eliminate_list_of_n_d1_d2(self, d1, d2, n_vals, bound_t=150):
        problematic_n = []
        for n in n_vals:
            success_n = self._eliminate_n_newforms_d1_d2(n, d1, d2, bound_t=bound_t)
            if not success_n:
                problematic_n.append(n)
        return problematic_n

    def _eliminate_n_newforms_d1_d2(self, n, d1, d2, bound_t=50):
        Ex = lambda x: EllipticCurve([0, 2 * x, 0, x ** 2 + 1, 0])
        pair_info = self.info[f"pair_{d1}_{d2}"]
        not_eliminated_newforms = pair_info['failed_newforms'].copy()
        for t in range(2, bound_t + 1):
            l = ZZ(t) * n + 1
            if (l in Primes()) and (ZZ(self.k / 2) % l != 0):
                for Ef in pair_info['failed_newforms']:
                    if Ef not in not_eliminated_newforms:
                        continue
                    suitable_Ef_l = True
                    Fl = FiniteField(l)
                    fkbar = self.fk.change_ring(Fl)
                    y = polygen(Fl, 'y')
                    t_unit_roots = [r[0] for r in (y ** t - 1).roots()]
                    alf = Ef.ap(l)
                    for zt in t_unit_roots:
                        x0s = [r[0] for r in (y ** 2 + 1 - (2 * Fl(d2) * zt) / Fl(d1)).roots()]
                        for x0 in x0s:
                            y2 = (fkbar(x0**2) * d2) / (Fl(2)**(self.k-2) * Fl(d1))
                            if y2 in t_unit_roots or y2.is_zero():
                                aEx0 = l + 1 - Ex(x0).order()
                                diff = (aEx0 - alf)
                                if diff % n == 0:
                                    suitable_Ef_l = False
                                    break
                        if not suitable_Ef_l:
                            break
                    if l % 4 == 1:
                        diff = (4 - alf ** 2)
                        if diff % n == 0:
                            suitable_Ef_l = False
                    if suitable_Ef_l:
                        not_eliminated_newforms.remove(Ef)
                if len(not_eliminated_newforms) == 0:
                    return True
        return False
