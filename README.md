## Sum of Consecutive Powers as a Perfect Power

This repository contains the code for the computations of the paper *'Sum of Consecutive Poweres as a Perfect Power'* 
written by Angelos Koutsianas and Nikolaos Tzanakis.


### Contents of the repository

- **elimination_step.py**: This file contains the class *'SumOfConsecutivePowersModularMethod'*. In this class the 
function *'eliminate_newforms_method_1'* is the implementation of the computations in Proposition 6.2. The function
*'eliminate_newforms_method_2'* is the implementation of the computations in Proposition 6.4.

- **linear_forms_in_logarithms.py**: This file contains the class *'LinearFormsInLogarithms'*. In this class the
function *'linear_forms_in_logarithms_best_bound'* which gives the upper bound n<sub>0</sub> according to Section 7 and 
Table 2.

- **small_exponents.py**: This file contains the class *'SmallExponentSolutions'*. In this class the function
*'solve_thue_equations'* solves, using Pari/GP, all the necessary Thue equation in Proposition 5.1.


### The computations

If you want to do/repeat the computations supporting the paper run the following in sage.

```commandline
import importlib
import time

cd /to_the_folder_of_the_repo/
ncpus = 10 # number of cpus your machine can support
```

#### Small exponents and Thue equations

The computations for Proposition 5.1.

```commandline
from small_exponents.small_exponents import SmallExponentSolutions

for k in [10, 26, 30, 34, 50, 58, 70, 74, 78, 82, 90]:
    print(f"#### k: {k} ####")
    for n in [3, 4, 5, 7]:
        ses = SmallExponentSolutions(ZZ(k), n)
        t0 = time.time()
        sols = ses.solve_thue_equations(ncpus=ncpus)
        t1 = time.time()
        print(f"n: {n}, time: {t1 - t0}, sols: {sols}")
```

#### Modular method

The computations for Propositions 6.2 and 6.4.

```commandline
from modular_method.elimination_step import SumOfConsecutivePowersModularMethod

bound_t = 1050
for k in [10, 26, 30, 34, 50, 58, 70, 74, 78, 82, 90]:
    print(f"k: {k}")
    t0 = time.time()
    scpmm = SumOfConsecutivePowersModularMethod(ZZ(k))
    scpmm.eliminate_newforms_method_1(primes_bound=50)
    t1 = time.time()
    print(f"Time for the eliminations step 1: {t1 - t0}.")
    scpmm.eliminate_newforms_method_2(bound_t=bound_t, ncpus=ncpus)
    print(f"")
```

#### Linear forms in logarithms

The computations for the upper bound n<sub>0</sub> of n in Table 2.

```commandline
from bound_n.linear_forms_in_logarithms import LinearFormsInLogarithms

for k in range(6, 100):
    if k % 4 == 2:
        lfil = LinearFormsInLogarithms(ZZ(k))
        n0 = lfil.linear_forms_in_logarithms_best_bound()
        print(f"k: {k}, n0: {n0}")
```



