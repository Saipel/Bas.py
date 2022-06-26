import numpy as np
import random

def Bas (Nsys):
    bas = np.empty(4 * Nsys, dtype=str)
    for i in range(0, 4*Nsys, 4):
        bas[0+i] = str("u")
        bas[1+i] = str(i/4+1)
        bas[2+i] = str("d")
        bas[3+i] = str(i/4+1)
    return (bas)

def sHam(bas, Nsys):
    ham = np.empty((2*Nsys, 2*4*Nsys), dtype=str)
    for i in range (2*Nsys):
        for j in range(0, 2*4*Nsys, 4):
            ham[i][j + 0] = str(bas[0 + 2*i])
            ham[i][j + 1] = str(bas[1 + 2*i])
            ham[i][j + 2] = str(bas[0 + 2*int(j/4)])
            ham[i][j + 3] = str(bas[1 + 2*int(j/4)])
    return(ham)

def out (ham, Nsys):
    for i in range (2*Nsys):
        for j in range(0, 2*4*Nsys, 4):
            print("<", ham[i][j], ham[i][j+1],"|", "H", "|", ham[i][j+2], ham[i][j+3], ">", end="\t")
        print("\n")

def zHam (ham, Nsys, w, t):
    fham = np.empty((2 * Nsys, 2 * Nsys), dtype=float)
    for i in range (2*Nsys):
        for j in range(0, 2*4*Nsys, 4):
            if ham[i][j+0] == ham[i][j+2] and ham[i][j+1] == ham[i][j+3]:
                fham[i][int(j/4)] = random.uniform(-w/2, w/2)

            if ham[i][j+0] == ham[i][j+2] and ham[i][j+1] == ham[i][j+3]:
                fham[i][int(j/4)] = random.uniform(-w/2, w/2)

            if ham[i][j+0] == ham[i][j+2] and ham[i][j+1] == ham[i][j+3]:
                fham[i][int(j/4)] = random.uniform(-w/2, w/2)

            if ham[i][j+0] == ham[i][j+2] and ham[i][j+1] == ham[i][j+3]:
                fham[i][int(j/4)] = random.uniform(-w/2, w/2)
    return(fham)