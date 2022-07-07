import threading as thg

import numpy as np

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
for i in range(Nlvl - 1):
    t[i] = w / 3

eps = np.empty(Nlvl)

Nsys = int(input("Nunber of sys\t"))

lamName = "Lam"
lam = []
Ipr = []
for i in range(-int((Nsys - 1) / 2), int(Nsys / 2) + 1):
    fham = Bas.zHam(ham, Nlvl, w, t, eps)  # заполнение массива

    # print(fham, end="\n\n")

    wb, vb = np.linalg.eigh(fham)
    # print(wb)
    wb2 = list(set(wb))
    del wb
    for i in range(len(wb2)):
        lam.append(wb2[i])

    for j in range(Nlvl):
        Ipr.append(Calc.IPR(vb, Nlvl, j))
    #print(Ipr)
    del vb
    # print(wb2)
    Output.OutPut(wb2, lamName)
    del fham
del ham

print("lam = ",lam)
print("Ipr = ",Ipr)

Calc.Sort(lam, Ipr, 0, len(lam)-1)

lam.append(0)
Ipr.append(0)

print("lam = ",lam)
print("Ipr = ",Ipr)

dosName = "DOS"
Ds = thg.Thread(target=Calc.DOS, args=(lamName, dosName, 4, 4, 0.0015, 0.001, Nsys))
Ds.start()
sobsName = "Sobs"
Calc.Energ(sobsName, 4, 4, 0.001)

aLam, aIpr = Calc.avengIpr(lam, Ipr)
print("aLam = ", aLam)
print("aIpr = ", aIpr)

avenLamName = "aLam"
avenIprName = "aIpr"

Output.OutPut(aLam, avenLamName)
Output.OutPut(aIpr, avenIprName)
