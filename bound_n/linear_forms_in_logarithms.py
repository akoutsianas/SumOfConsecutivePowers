from sage.all import log, polygen, prod, Subsets, Infinity, QQ, ZZ, RR, round, previous_prime, Primes

from config import laurent_constants


class LinearFormsInLogarithms:

    def __init__(self, k):
        if k % 4 != 2:
            raise ValueError(f"k is not equal to 2 mod(4)!")
        self.k = k
        self._fk = self._compute_fk()
        self._hks = self._compute_hks()
        self._hk = self._hks[self.k/2-2]
        self._lower_bound = 2
        self.x0 = self._compute_x0()

    def _compute_fk(self):
        x = polygen(QQ, 'x')
        g = ((x - 1)**self.k + (x + 1)**self.k)/2
        f = x.parent(g/(x**2 + 1))
        f = sum([fi*x**i for i, fi in enumerate(f.coefficients())])
        return f

    def _compute_hks(self):
        x = self._fk.parent().gen()
        h = self._fk(x - 1)
        hks = h.coefficients()
        return hks

    def _compute_x0(self):
        x0 = sum([abs(hi) for hi in self._hks[:-2]]) / self._hk
        return x0

    def _compute_ds(self):
        ds = []
        for d in ZZ(self.k/2).divisors():
            if len([p for p in ZZ(d).prime_factors() if p % 4 == 3]) == 0 and d != 1:
                ds.append(d)
        return ds

    def _compute_A_B(self, d):
        """
            We assume that tp = 1 for all prime divisors of d
        """
        S = set(d.prime_factors())
        pairs = []
        for P1 in Subsets(S):
            P2 = S.difference(P1)
            A = prod([p**(self.k * d.valuation(p)/2) for p in P1]) / prod([p**(self.k * d.valuation(p)/2) for p in P2])
            B = prod([p for p in P2]) / prod([p**(self.k/2 - 2) for p in P1])
            pairs.append((A, B))
        return pairs

    def linear_forms_in_logarithms_best_bound(self, step_n=1000):
        bound = {}
        ds = self._compute_ds()
        for d in ds:
            bound_d = Infinity
            for key, value in laurent_constants.items():
                b0 = self._linear_forms_in_logarithms_bound(d, key, value['c1'], step_n=step_n)
                bound_d = min([bound_d, b0])
            bound[d] = bound_d
        return bound

    def _linear_forms_in_logarithms_bound(self, d, m, C1, step_n=1000, y1_lbound=3):
        """
        INPUT:
            m: the constant in Laurent's paper in Table 1
            C1: the constant in Laurent's paper in Table 1

        OUTPUT:
            An upper bound of n using linear forms in logarithms.
        """

        pairs = self._compute_A_B(d)
        bound = 5
        for A, B in pairs:
            a2 = A/2**((self.k/2) - 2)
            logA2 = max(log(a2.height()), 1)
            C = log(1.0051 * self._hk / d)
            E = max([1/log(y1_lbound), self.k/2 - 1, log(B)/log(y1_lbound) + log(1.5)/log(y1_lbound) + self.k/2 - 1])
            bound_pair = 5

            ### Max is m ###
            temp_bound = round(RR(C1 * m**2 * logA2 * E + C/log(y1_lbound)))
            bound_pair = max([temp_bound, bound_pair])

            ### Max (log(b' + 0.21))
            logf = lambda x: C1 * E * logA2 * (log(x/logA2 + 1) + 0.21)**2 + C/log(y1_lbound)
            n = 5
            while n <= logf(n):
                n += step_n
            while n > logf(n):
                n -= 1
            bound_pair = ZZ(round(max([n, bound_pair])))
            bound = max([bound, bound_pair])
        if not bound in Primes():
            bound = previous_prime(bound)
        return bound

