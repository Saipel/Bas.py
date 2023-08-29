import numpy as np
#from numba import njit, prange

import Calc_with_inter as CWI

# основные константы
coulomb_potential = [4,4]
zone_wight = int(input("zone_wight\t"))
#number_of_case = int(input("number_of_case\t"))
systems_of_number = int(input("systems_of_number\t"))


# создание массива для записи значения энергии состояния и коэффициентов перед дельта функцией
array_with_energy_for_calculate = np.zeros((3, systems_of_number * 18))

number_of_elements_per_line = 0

for system_number in range(1, systems_of_number + 1):

    # создание параметров гамильтониана
    site_potential, jump_parameter = CWI.create_hamiltonian_parametrs(zone_wight)
   # print("site_potential = ", site_potential)
   # print("jump_parameter = ", jump_parameter)


    # создание нулевых для всех ground state
    one_part_ham = CWI.create_null_array_for_all_cases(1)
    two_part_ham = CWI.create_null_array_for_all_cases(2)
    three_part_ham = CWI.create_null_array_for_all_cases(3)
    four_part_ham = CWI.create_null_array_for_all_cases(4)


    # Запись параметров гамильтониана в нулевой массив
    one_part_ham = CWI.hamiltanian_arays_with_int_for_all_cases(1, one_part_ham, site_potential, jump_parameter,
                                                 coulomb_potential)
    two_part_ham = CWI.hamiltanian_arays_with_int_for_all_cases(2, two_part_ham, site_potential, jump_parameter,
                                                 coulomb_potential)
    #print("Hamiltonian: \n", two_part_ham)
    three_part_ham = CWI.hamiltanian_arays_with_int_for_all_cases(3, three_part_ham, site_potential, jump_parameter,
                                                 coulomb_potential)
    four_part_ham = CWI.hamiltanian_arays_with_int_for_all_cases(4, four_part_ham, site_potential, jump_parameter,
                                                 coulomb_potential)



    # вычисление массивов собственных значений для всех случаев
    one_eigenvalue_array = CWI.calc_eigenvalue_array(one_part_ham, 1)
    two_eigenvalue_array = CWI.calc_eigenvalue_array(two_part_ham, 2)
   # print("Eigenvalue: ", two_eigenvalue_array)
    three_eigenvalue_array = CWI.calc_eigenvalue_array(three_part_ham, 3)
    four_eigenvalue_array = CWI.calc_eigenvalue_array(four_part_ham, 4)

    # вычисление собственных векторов
    one_eigenvectors_array = CWI.calc_eigenvectors_of_an_array(one_part_ham)
    two_eigenvectors_array = CWI.calc_eigenvectors_of_an_array(two_part_ham)
    #print("Eigenvector: \n", two_eigenvectors_array)
    three_eigenvectors_array = CWI.calc_eigenvectors_of_an_array(three_part_ham)

    # Составление массива энергий для дельта функций
    array_with_state_energy = CWI.array_with_state_energy(one_eigenvalue_array, two_eigenvalue_array, three_eigenvalue_array, four_eigenvalue_array)


    # комментарий
    number_of_elements_per_line = CWI.array_with_self_energ_and_(array_with_state_energy, one_eigenvectors_array, two_eigenvectors_array, three_eigenvectors_array,
                                   array_with_energy_for_calculate, number_of_elements_per_line, system_number)


CWI.quick_sort(array_with_energy_for_calculate, 0, len(array_with_energy_for_calculate[0]) - 1)
print(CWI.binary_search_recursive(array_with_energy_for_calculate, 9.32, 0,len(array_with_energy_for_calculate[0]) - 1))

#CWI.array_out(array_with_energy_for_calculate)

# just a hyper threading
import threading as thg

x = thg.Thread(target=CWI.DOS_calc_for_inter_system, args=(array_with_energy_for_calculate, zone_wight, systems_of_number,coulomb_potential))
y = thg.Thread(target=CWI.ensemble_averaged_GIPR, args=(array_with_energy_for_calculate, zone_wight, systems_of_number,coulomb_potential))
z = thg.Thread(target=CWI.energy_distribution, args=(zone_wight, array_with_energy_for_calculate, systems_of_number,coulomb_potential))
x.start()
y.start()
z.start()
