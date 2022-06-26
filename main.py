import numpy as np
import cmath
import Bas

Nsys = int(input("Num of sys\t"))
bas = Bas.Bas(Nsys)
print(bas)
ham = Bas.sHam(bas, Nsys)
Bas.out(ham,Nsys)

w = r