from itertools import product
from sage.all import polygen, QuadraticField, prod, gp, ZZ
from sage.parallel.use_fork import p_iter_fork


class SmallExponentSolutions:

    def __init__(self, k, n):
        self.k = k
        self.n = n
        self._K = QuadraticField(-1, 'i')
        self._ri = self._K.gen()
        self._conj = [emb for emb in self._K.embeddings(self._K) if emb(self._ri) != self._ri][0]
        self._alphas = self._compute_alphas()
        self._us = self._compute_us()

    def _compute_alphas(self):
        primes = [p for p in self.k.prime_factors() if p % 4 == 1]
        if len(primes) == 0:
            return [1]
        primes_gen = [self._K.prime_above(p).gens_reduced()[0] for p in primes]
        alphas = []
        for v in product(range(self.n), repeat=len(primes_gen)):
            a = prod([ai**ei for ai, ei in zip(primes_gen, v)])
            alphas.append(a)
        alphas += [self._conj(a) for a in alphas if self._conj(a) != a]
        return alphas

    def _compute_us(self):
        if self.n != 4:
            return [1 + self._ri]
        else:
            return [1 + self._ri, -1 + self._ri, 1 - self._ri, -1 - self._ri]

    def solve_thue_equations(self, ncpus=1):
        betas = [(1 + self._ri) * alpha for alpha in self._alphas]
        fork_iterator = p_iter_fork(ncpus=ncpus)
        inputs = [
            (
                [[beta for j, beta in enumerate(betas) if j % ncpus == i]],
            ) for i in range(ncpus)]
        sols_thue = list(fork_iterator(self._solve_thue_equations_for_pairs, inputs))
        sols = []
        for sol in sols_thue:
            for x0 in sol[1]:
                if x0 not in sols:
                    sols.append(x0)
        return sols

    def _solve_thue_equations_for_pairs(self, betas):
        sols = []
        x = polygen(self._K, 'x')
        for beta in betas:
            fr = (beta * (x + self._ri) ** self.n + self._conj(beta) * (x - self._ri) ** self.n) / 2
            fc = (beta * (x + self._ri) ** self.n - self._conj(beta) * (x - self._ri) ** self.n) / (2 * self._ri)
            fr = fr.change_ring(ZZ)
            fc = fc.change_ring(ZZ)
            sols_thue = gp.thue(gp.thueinit(fr, 1), 1).sage()
            sols_thue += gp.thue(gp.thueinit(fr, 1), -1).sage()
            sols_thue += gp.thue(gp.thueinit(fc, 1), 1).sage()
            sols_thue += gp.thue(gp.thueinit(fc, 1), -1).sage()
            for t0, s0 in sols_thue:
                x0 = (beta * (t0 + s0 * self._ri) ** self.n + self._conj(beta) * (t0 - s0 * self._ri) ** self.n) / 2
                x0 = ZZ(x0)
                w = ((x0 - 1) / 2) ** self.k + ((x0 + 1) / 2) ** self.k
                if w.is_perfect_power():
                    if x0 not in sols:
                        sols.append(x0)
        return sols

