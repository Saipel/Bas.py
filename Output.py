def OutPut(mas, fname):
    fail = open(fname + ".txt", 'a')  # Открытие файла для вывода
    for i in range(len(mas)):  # Цикл для оследовательного вывода обственных значений
        fail.write(str(mas[i]) + "\n")
    #fail.write()
    fail.close()