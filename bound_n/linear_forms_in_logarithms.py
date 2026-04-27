from sage.all import (log, polygen, prod, Subsets, Infinity, QQ, ZZ, RR, round, previous_prime, Primes, exp, gcd, sqrt,
                      primes_first_n)

from bound_n.laurent_parameters import laurent_constants


class LinearFormsInLogarithms:

    def __init__(self, k):
        if k % 4 != 2:
            raise ValueError(f"k is not equal to 2 mod(4)!")
        self.k = k
        self._fk = self._compute_fk()
        self._hks = self._compute_hks()
        self._hkm2 = self._hks[self.k/2-2]
        self._lower_bound = 2
        self.x0 = self._compute_x0()
        self._c1 = self._compute_c1()
        self._A_upper = self.k**(self.k/2) / 2
        self._n1 = self._compute_n1()
        self._n2 = self._compute_n2()
        self._n3 = self._compute_n3()

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
        x0 = int(max([sum([abs(hi) for hi in self._hks[:-2]]) / self._hkm2, 200*self._hkm2, self.k])) + 1
        return x0

    def _solutions_for_y1_one(self):
        for d2 in ZZ(self.k/2).divisors():
            x = 2 * d2 - 1
            if x.is_square() and x != 1:
                x = ZZ(sqrt(x))
                yn = ZZ(((x - 1)/2)**self.k + ((x + 1)/2)**self.k)
                if yn.is_perfect_power():
                    print(f"We have a solution for x={x} and y1=1.")

    def _compute_c1(self):
        c1 = 0.50255 * self.k**2 * (self.k/2 - 1)
        return c1

    def _compute_n1(self):
        n1 = int(log(self.k * self.x0/4) / log(3))
        return n1

    def _compute_n2(self):
        n2 = max([log(self._c1/log(self._A_upper))/log(3), log(self._c1)/log(3), (log(self._A_upper) + 1)/log(1.5)])
        n2 = int(n2)
        return n2

    def _compute_n3(self):
        n3 = int(100 * (self.k * log(self.k)) / (log(1.5 * 3**(self.k/2 - 1))))
        return n3

    def linear_forms_in_logarithms_best_bound(self, step_n=1000):
        n0 = min([self._linear_forms_in_logarithms(m, C2) for m, C2 in laurent_constants.items()])
        return n0

    def _linear_forms_in_logarithms(self, m, C2, step_n=1000):
        n4 = self._compute_n4(m)
        w1 = (self.k/2 - 1 + log(1.5)/log(3)) * self.k * log(self.k)/2
        w2 = log(0.50255 * self.k**2 * (self.k/2 - 1)) / log(3)
        f = lambda n: C2 * (log((2.01*n)/(self.k*log(self.k))))**2 * w1 + w2 - n

        n0 = 7
        while f(n0) > 0:
            n0 += step_n
        while f(n0) <= 0:
            n0 -= 1
        n0 = int(max([self._n1, self._n2, self._n3, n4, n0]))
        return n0

    def _compute_n4(self, m):
        n4 = int(exp(m - 0.38) * self.k * log(self.k) / 2)
        return n4

