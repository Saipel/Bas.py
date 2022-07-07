import cmath
import numpy as np

# Функция расчета DOS и вывода значений в файл
# Сравнить с существующей программой
def DOS(name, dname, min=None, max=None, gam=None, step=None, Nsys=None):
    dfail = open(dname + ".txt", "w")

    gr = - float(min)
    while gr < float(max):
        fail = open(name + ".txt", "r")
        dos = 0
        for line in fail:
            dos = dos + (1 / (cmath.pi * Nsys)) * (gam / (gam ** 2 + (float(line) - gr) ** 2))
        dfail.write(str(dos) + "\n")
        gr += step
        fail.close()

    dfail.close()


# Функция вывода значений энергии для построения графиков DOS(E)
# Работает корректно
def Energ(dname, min=None, max=None, step=None):
    dfail = open(dname + ".txt", "w")
    gr = - float(min)
    while gr < float(max):
        dfail.write(str(gr) + "\n")
        gr += step
    dfail.close()


# Функция расчета IPR для одной системы
# Необходимо проверить корректность работы (сравнить с существующей программой)
def IPR(mas, Nlvl, n):
    ipr = 0
    nn = 0  # Вспомогательная переменная

    for i in range(Nlvl):
        ipr += mas[n * 2][i] ** 2 + mas[n * 2 + 1][i] ** 2
        nn += (mas[n * 2][i] + mas[n * 2 + 1][i]) ** 2

    ipr = ipr / nn
    return ipr


# Функция быстрой сортировки двух массивов
# Работает корректно
def Sort(mas1, mas2, fst, lst):
    if fst >= lst:
        return mas1, mas2

    i, j = fst, lst
    pivot = mas1[int((fst + lst) / 2)]
    while i <= j:
        while mas1[i] < pivot:
            i += 1
        while mas1[j] > pivot:
            j -= 1
        if i <= j:
            mas1[i], mas1[j] = mas1[j], mas1[i]
            mas2[i], mas2[j] = mas2[j], mas2[i]
            i, j = i + 1, j - 1

    else:
        return Sort(mas1, mas2, fst, j), Sort(mas1, mas2, i, lst)


def avengIpr(Lam, Ipr):
    j = 0
    sch = 0

    aLam = np.zeros(len(Lam) - 1)      #Avenger lam
    aIpr = np.zeros(len(Ipr) - 1)      #Avenger IPR

    i = 0
    while i < len(Lam) - 1:
        j = 0
        if int(Lam[i + j] * 1000) == int(Lam[i + j + 1] * 1000):
            while int(Lam[i + j] * 1000) == int(Lam[i + j + 1] * 1000):
                aLam[sch] = Lam[i + j]
                aIpr[sch] = aIpr[sch] + Ipr[i + j]
                j += 1
            aIpr[sch] = aIpr[sch]+Ipr[i+j]
            aIpr[sch] = aIpr[sch]/(j+1)

        else:
            aLam[sch] = Lam[i]
            aIpr[sch] = Ipr[i]

        sch += 1
        i = i + j + 1

    return aLam, aIpr
    del aLam
    del aIpr
