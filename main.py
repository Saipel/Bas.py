import threading as thg

import numpy as np

from numba import prange, jit

import Bas
import Calc
import Output

Nlvl = int(input("Num of lvl\t"))
bas = Bas.Bas(Nlvl)

print(bas)

ham = Bas.sHam(bas, Nlvl)
del bas
Bas.out(ham, Nlvl)

w = int(input("Zone width\t"))
t = np.empty(Nlvl - 1)

# Рассматривается только перескоки с соседних уровней
for i in prange(Nlvl - 1):
    t[i] = w / 6

eps = np.empty(Nlvl)

Nsys = int(input("Number of sys\t"))

lamName = "Lam"
lam = []
Ipr = []
for i in prange(Nsys):
    fham = Bas.zHam(ham, Nlvl, w, t, eps)  # заполнение массива

    # print(fham, end="\n\n")

    wb, vb = np.linalg.eigh(fham)
    #print(wb)
    k = 0
    while k < len(wb):
        lam.append(wb[k])
        k+=2
    del wb
    #print(wb2)

    for j in range(Nlvl):
        Ipr.append(Calc.IPR(vb, Nlvl, j))
    #print(Ipr)
    del vb
    # print(wb2)

    del fham
del ham

Output.OutPut(lam, lamName)

#print("lam = ",lam)
#print("Ipr = ",Ipr)

Calc.Sort(lam, Ipr, 0, len(lam)-1)

lam.append(0)
Ipr.append(0)



#print("lam = ",lam)
#print("Ipr = ",Ipr)

dosName = "DOS"
Ds = thg.Thread(target=Calc.DOS, args=(lamName, dosName, lam, 4, 4, 0.0015, 0.001, Nsys))
Ds.start()
sobsName = "Sobs"
Calc.Energ(sobsName, 4, 4, 0.001)

aLam, aIpr = Calc.avengIpr(lam, Ipr)

###################################
aLam, aIpr = Calc.Optim(aLam,aIpr)

avenLamName = "aLam"
avenIprName = "aIpr"

x = thg.Thread(target=Output.OutPut, args = (aLam, avenLamName))
y = thg.Thread(target=Output.OutPut, args = (aIpr, avenIprName))
x.start()
y.start()