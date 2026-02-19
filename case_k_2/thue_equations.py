K.<ri> = QuadraticField(-1)
a = polygen(K, 'a')

def thue_equation_p(p):
    f = ((1 - ri) * (a - ri)**p + (1 + ri) * (a + ri)**p)/2
    f = f.change_ring(ZZ)
    tnf = gp.thueinit(f, 1)
    sols = gp.thue(tnf, 1).sage()
    return sols


exceptional_primes = []
for p in prime_range(3, 500):
    sols = thue_equation_p(p)
    for sol in sols:
        y = sol[0]**2 + sol[1]**2
        if y != 1:
            print(f"We have got a solution different to y=1 for p={p}")
            exceptional_primes.append(p)
    print(f"The only solutions for p={p} is for y=1!")
print(f"The exceptional primes are {exceptional_primes}")
