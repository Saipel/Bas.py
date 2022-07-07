import numpy as np
import random

def Bas(Nsys):
    bas = np.empty(4 * Nsys, dtype=str)
    for i in range(0, int(4 * Nsys/2), 2):
        bas[i] = str("u")
        bas[i + 1] = str(int(i/2)+1)
    j = 0
    for i in range(int(4 * Nsys/2), 4*Nsys, 2):
        bas[i] = str("d")
        bas[i + 1] = str(j+1)
        j+=1
    return (bas)
    del bas



def sHam(bas, Nsys):
    ham = np.empty((2*Nsys, 2*4*Nsys), dtype=str)
    for i in range (2*Nsys):
        for j in range(0, 2*4*Nsys, 4):
            ham[i][j + 0] = str(bas[0 + 2*i])
            ham[i][j + 1] = str(bas[1 + 2*i])
            ham[i][j + 2] = str(bas[0 + 2*int(j/4)])
            ham[i][j + 3] = str(bas[1 + 2*int(j/4)])
    return(ham)
    del ham

def out (ham, Nsys):
    for i in range (2*Nsys):
        for j in range(0, 2*4*Nsys, 4):
            print("<", ham[i][j], ham[i][j+1],"|", "H", "|", ham[i][j+2], ham[i][j+3], ">", end="\t")
        print("\n")

def zHam (ham, Nsys, w, t,eps):

#    Заолняем массив потенцилов
    for i in range(len(eps)):
        eps[i] = random.uniform(-w/2, w/2)

    fham = np.zeros((2 * Nsys, 2 * Nsys), dtype=float)
    p = 0       #переменная отвечает за перебор значений в массиве потенциалов
    for i in range(2*Nsys):
        k = 0               #переменая отвечает за перебор элементов в массивве перескока

        if p == len(eps):
            p = 0

        for j in range(0, 2*4*Nsys, 4):
            if ham[i][j+0] == ham[i][j+2] and ham[i][j+1] == ham[i][j+3]:
                fham[i][int(j/4)] = eps[p]
                p += 1
            if ham[i][j+0] == ham[i][j+2] and (int(ham[i][j+1]) == (int(ham[i][j+3]) + 1)
                or int(ham[i][j+1]) + 1 == int(ham[i][j+3]) or int(ham[i][j+1]) == int(ham[i][j+3])+Nsys-1
                or int(ham[i][j+1]) + Nsys- 1 == int(ham[i][j+3])):
                fham[i][int(j/4)] = -t[k]
                k += 1
    #print(fham)
    return(fham)
    del fham