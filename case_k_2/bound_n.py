f = lambda p: 2*log(1.004*sqrt(2)) / log(5) + 72.8 * log(p+1)**2 - p

p = 7
while f(p) >= 0:
    p += 1

p -= 1
print(f"It should hold p<= {p}")