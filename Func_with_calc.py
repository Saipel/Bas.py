import threading as thg

import numpy as np

from numba import prange, jit

import Bas
import Calc
import Output

def two_state_sys_with_inter():
    return 0

def sys_without_inter_for_diff_number_state(ham, Nlvl, w, t, array_of_potential, lam, Ipr):
    fham = Bas.zHam(ham, Nlvl, w, t, array_of_potential)  # заполнение массива

    # print(fham, end="\n\n")

    wb, vb = np.linalg.eigh(fham)  # Calc of eigenvalues and eigenvectors
    # print(wb)

    # Create eigenvalue array for use
    k = 0
    while k < len(wb):
        lam.append(wb[k])
        k += 2
    del wb
    # print(wb2)

    # Create IPR array for next calc
    for j in range(Nlvl):
        Ipr.append(Calc.IPR(vb, Nlvl, j))
    # print(Ipr)
    del vb
    # print(wb2)

    del fham

    return 0

def all_output_for_sys_without_inter(lam, lamName, Ipr, Nsys):
    Output.OutPut(lam, lamName)

    # print("lam = ",lam)
    # print("Ipr = ",Ipr)

    Calc.Sort(lam, Ipr, 0, len(lam) - 1)

    lam.append(0)
    Ipr.append(0)

    # print("lam = ",lam)
    # print("Ipr = ",Ipr)

    dosName = "DOS"
    Ds = thg.Thread(target=Calc.DOS, args=(lamName, dosName, lam, 4, 4, 0.0015, 0.001, Nsys))
    Ds.start()
    sobsName = "Sobs"
    Calc.Energ(sobsName, 4, 4, 0.001)

    aLam, aIpr = Calc.avengIpr(lam, Ipr)

    ###################################
    aLam, aIpr = Calc.Optim(aLam, aIpr)

    avenLamName = "aLam"
    avenIprName = "aIpr"

    x = thg.Thread(target=Output.OutPut, args=(aLam, avenLamName))
    y = thg.Thread(target=Output.OutPut, args=(aIpr, avenIprName))
    x.start()
    y.start()