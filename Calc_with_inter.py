import math
import random

import numpy as np


# from numba import njit, prange


def create_null_array_for_all_cases(number_of_case):
    if number_of_case == 1:
        null_array = np.zeros((4, 4))
    if number_of_case == 2:
        null_array = np.zeros((6, 6))
    if number_of_case == 3:
        null_array = np.zeros((4, 4))
    if number_of_case == 4:
        null_array = 0

    return null_array


def create_hamiltonian_parametrs(zone_wight):
    site_potential = np.zeros(2)
    for i in range(len(site_potential)):
        site_potential[i] = random.uniform(-zone_wight / 2, zone_wight / 2)

    jump_parameter = zone_wight / 6

    return site_potential, jump_parameter


def hamiltanian_arays_with_int_for_all_cases(number_of_case, null_array, site_potential, jump_parameter,
                                             coulomb_potential):
    if number_of_case == 1:
        null_array[0][0] = site_potential[0]
        null_array[1][1] = site_potential[1]
        null_array[2][2] = site_potential[0]
        null_array[3][3] = site_potential[1]
        null_array[0][1] = null_array[1][0] = -jump_parameter
        null_array[2][3] = null_array[3][2] = -jump_parameter

    if number_of_case == 2:
        null_array[0][0] = null_array[1][1] = null_array[2][2] = null_array[3][3] = site_potential[0] + site_potential[
            1]
        null_array[4][4] = 2 * site_potential[0] + coulomb_potential
        null_array[5][5] = 2 * site_potential[1] + coulomb_potential
        null_array[3][4] = null_array[3][5] = null_array[4][3] = null_array[5][3] = - math.sqrt(2) * jump_parameter

    if number_of_case == 3:
        null_array[0][0] = 2 * site_potential[0] + site_potential[1]
        null_array[1][1] = 2 * site_potential[1] + site_potential[0]
        null_array[2][2] = 2 * site_potential[0] + site_potential[1]
        null_array[3][3] = 2 * site_potential[1] + site_potential[0]
        null_array[0][1] = null_array[1][0] = -jump_parameter
        null_array[2][3] = null_array[3][2] = -jump_parameter

    if number_of_case == 4:
        null_array = 2 * site_potential[0] + 2 * site_potential[1]

    return null_array


# Need to finish of JIT
# @njit()
def approx_of_delta_function(energ_resolution, energ_of_state, energy_distribution):
    approx_value = 1 / math.pi * (
            energ_resolution / ((energ_of_state - energy_distribution) ** 2 + energ_resolution ** 2))

    return approx_value


def calc_eigenvalue_array(array_for_calc, number_of_case):
    array_for_output = []
    if number_of_case != 4:
        intermediate_stage, _ = np.linalg.eig(array_for_calc)
        # print("intermediate_stage = ", intermediate_stage)
        return intermediate_stage
    else:
        array_for_output.append(array_for_calc)

    return array_for_output


def calc_eigenvectors_of_an_array(array_for_calc):
    _, intermediate_stage = np.linalg.eig(array_for_calc)
    # print("intarmidate of vector array = \n", intermediate_stage)

    return intermediate_stage


# Расчет массива энергий состояний
def array_with_state_energy(one_part_array, two_part_array, three_part_array, four_part_array):
    intermediate_stage = np.zeros(18)

    # zero particle ground state
    intermediate_stage[0] = one_part_array[0]
    intermediate_stage[1] = one_part_array[1]

    # one particle ground state
    intermediate_stage[2] = two_part_array[5] - one_part_array[0]
    intermediate_stage[3] = one_part_array[0]
    intermediate_stage[4] = two_part_array[0] - one_part_array[0]
    intermediate_stage[5] = two_part_array[1] - one_part_array[0]
    intermediate_stage[6] = two_part_array[2] - one_part_array[0]

    # two particle ground state
    intermediate_stage[7] = three_part_array[0] - two_part_array[0]
    intermediate_stage[8] = three_part_array[1] - two_part_array[0]
    intermediate_stage[9] = two_part_array[0] - one_part_array[0]
    intermediate_stage[10] = two_part_array[0] - one_part_array[1]

    # three particle ground state
    intermediate_stage[11] = three_part_array[0] - two_part_array[5]
    intermediate_stage[12] = four_part_array - three_part_array[0]
    intermediate_stage[13] = three_part_array[0] - two_part_array[0]
    intermediate_stage[14] = three_part_array[0] - two_part_array[1]
    intermediate_stage[15] = three_part_array[0] - two_part_array[2]

    # four particle ground state
    intermediate_stage[16] = four_part_array - three_part_array[0]
    intermediate_stage[17] = four_part_array - three_part_array[1]

    return intermediate_stage


# Допили название быдло...
# @njit(parallel=True)
def array_with_self_energ_and_(energ_of_state, one_part_ham, two_part_ham,
                               three_part_ham, array_with_energy_for_calculate, number_of_elements_per_line,
                               syst_number):
    for i in range(2):
        array_with_energy_for_calculate[0][i + number_of_elements_per_line] = energ_of_state[
            i + number_of_elements_per_line - 18 * (syst_number - 1)]

    # energy first translation on diff site
    array_with_energy_for_calculate[1][0 + number_of_elements_per_line] = 0.5 * (
                one_part_ham[0][0] + one_part_ham[2][0]) ** 2
    array_with_energy_for_calculate[2][0 + number_of_elements_per_line] = 0.5 * (
                one_part_ham[1][0] + one_part_ham[3][0]) ** 2

    # energy first translation on diff site
    array_with_energy_for_calculate[1][1 + number_of_elements_per_line] = 0.5 * (
                one_part_ham[0][2] + one_part_ham[2][2]) ** 2
    array_with_energy_for_calculate[2][1 + number_of_elements_per_line] = 0.5 * (
                one_part_ham[1][2] + one_part_ham[3][2]) ** 2

    number_of_elements_per_line += 2

    # Look in article
    for i in range(5):
        array_with_energy_for_calculate[0][i + number_of_elements_per_line] = energ_of_state[
            i + number_of_elements_per_line - 18 * (syst_number - 1)]

    # energy first translation on diff site
    array_with_energy_for_calculate[1][0 + number_of_elements_per_line] = 0.25 * 3.0 / 2.0 * (
                one_part_ham[0][0] + one_part_ham[2][0]) ** 2
    array_with_energy_for_calculate[2][0 + number_of_elements_per_line] = 0.25 * 3.0 / 2.0 * (
                one_part_ham[1][0] + one_part_ham[3][0]) ** 2

    # energy second translation on diff site
    array_with_energy_for_calculate[1][1 + number_of_elements_per_line] = 0.25 * (
                one_part_ham[0][0] + one_part_ham[2][0]) ** 2
    array_with_energy_for_calculate[2][1 + number_of_elements_per_line] = 0.25 * (
                one_part_ham[1][0] + one_part_ham[3][0]) ** 2

    # energy second translation on diff site
    array_with_energy_for_calculate[1][2 + number_of_elements_per_line] = 0.25 * (
            1 / math.sqrt(2) * two_part_ham[3][0] * (one_part_ham[1][0] + one_part_ham[3][0]) +
            two_part_ham[4][0] * (one_part_ham[0][0] + one_part_ham[2][0])) ** 2

    array_with_energy_for_calculate[2][2 + number_of_elements_per_line] = 0.25 * (
            1 / math.sqrt(2) * two_part_ham[3][0] * (one_part_ham[0][0] + one_part_ham[2][0]) +
            two_part_ham[5][0] * (one_part_ham[1][0] + one_part_ham[3][0])) ** 2

    # energy second translation on diff site
    array_with_energy_for_calculate[1][3 + number_of_elements_per_line] = 0.25 * (
            1 / math.sqrt(2) * two_part_ham[3][2] * (one_part_ham[1][0] + one_part_ham[3][0]) +
            two_part_ham[4][2] * (one_part_ham[0][0] + one_part_ham[2][0])) ** 2

    array_with_energy_for_calculate[2][3 + number_of_elements_per_line] = 0.25 * (
            1 / math.sqrt(2) * two_part_ham[3][2] * (one_part_ham[0][0] + one_part_ham[2][0]) +
            two_part_ham[5][2] * (one_part_ham[1][0] + one_part_ham[3][0])) ** 2

    # energy second translation on diff site
    array_with_energy_for_calculate[1][4 + number_of_elements_per_line] = 0.25 * (
            1 / math.sqrt(2) * two_part_ham[3][1] * (one_part_ham[1][0] + one_part_ham[3][0]) +
            two_part_ham[4][1] * (one_part_ham[0][0] + one_part_ham[2][0])) ** 2

    array_with_energy_for_calculate[2][4 + number_of_elements_per_line] = 0.25 * (
            1 / math.sqrt(2) * two_part_ham[3][1] * (one_part_ham[0][0] + one_part_ham[2][0]) +
            two_part_ham[5][1] * (one_part_ham[1][0] + one_part_ham[3][0])) ** 2

    number_of_elements_per_line += 5

    for i in range(4):
        array_with_energy_for_calculate[0][i + number_of_elements_per_line] = energ_of_state[
            i + number_of_elements_per_line - 18 * (syst_number - 1)]

        # energy first translation on diff site
    array_with_energy_for_calculate[1][0 + number_of_elements_per_line] = 0.5 * (
            - 1 / math.sqrt(2) * two_part_ham[3][0] * (three_part_ham[0][0] + three_part_ham[2][0]) +
            two_part_ham[5][0] * (three_part_ham[1][0] + three_part_ham[3][0])) ** 2

    array_with_energy_for_calculate[2][0 + number_of_elements_per_line] = 0.5 * (
            - 1 / math.sqrt(2) * two_part_ham[3][0] * (three_part_ham[1][0] + three_part_ham[3][0]) +
            two_part_ham[4][0] * (three_part_ham[0][0] + three_part_ham[2][0])) ** 2

    # energy first translation on diff site
    array_with_energy_for_calculate[1][1 + number_of_elements_per_line] = 0.5 * (
            - 1 / math.sqrt(2) * two_part_ham[3][0] * (three_part_ham[0][2] + three_part_ham[2][2]) +
            two_part_ham[5][0] * (three_part_ham[1][2] + three_part_ham[3][2])) ** 2

    array_with_energy_for_calculate[2][1 + number_of_elements_per_line] = 0.5 * (
            - 1 / math.sqrt(2) * two_part_ham[3][0] * (one_part_ham[1][2] + one_part_ham[3][2]) +
            two_part_ham[4][0] * (one_part_ham[0][2] + one_part_ham[2][2])) ** 2

    # energy first translation on diff site
    array_with_energy_for_calculate[1][2 + number_of_elements_per_line] = 0.5 * (
            1 / math.sqrt(2) * two_part_ham[3][0] * (one_part_ham[0][0] + one_part_ham[2][0]) +
            two_part_ham[5][0] * (one_part_ham[1][0] + one_part_ham[3][0])) ** 2

    array_with_energy_for_calculate[2][2 + number_of_elements_per_line] = 0.5 * (
            1 / math.sqrt(2) * two_part_ham[3][0] * (one_part_ham[1][0] + one_part_ham[3][0]) +
            two_part_ham[4][0] * (one_part_ham[0][0] + one_part_ham[2][0])) ** 2

    # energy first translation on diff site
    array_with_energy_for_calculate[1][3 + number_of_elements_per_line] = 0.5 * (
            1 / math.sqrt(2) * two_part_ham[3][0] * (one_part_ham[0][2] + one_part_ham[2][2]) +
            two_part_ham[5][0] * (one_part_ham[1][2] + one_part_ham[3][2])) ** 2

    array_with_energy_for_calculate[2][3 + number_of_elements_per_line] = 0.5 * (
            1 / math.sqrt(2) * two_part_ham[3][0] * (one_part_ham[1][2] + one_part_ham[3][2]) +
            two_part_ham[4][0] * (one_part_ham[0][2] + one_part_ham[2][2])) ** 2

    number_of_elements_per_line += 4

    # Look in article
    for i in range(5):
        array_with_energy_for_calculate[0][i + number_of_elements_per_line] = energ_of_state[
            i + number_of_elements_per_line - 18 * (syst_number - 1)]

        # energy first translation on diff site
    array_with_energy_for_calculate[1][0 + number_of_elements_per_line] = 0.25 * 3.0 / 2.0 * (
            three_part_ham[0][0] + three_part_ham[2][0]) ** 2
    array_with_energy_for_calculate[2][0 + number_of_elements_per_line] = 0.25 * 3.0 / 2.0 * (
            three_part_ham[1][0] + three_part_ham[3][0]) ** 2

    # energy second translation on diff site
    array_with_energy_for_calculate[1][1 + number_of_elements_per_line] = 0.25 * (
            three_part_ham[0][0] + three_part_ham[2][0]) ** 2
    array_with_energy_for_calculate[2][1 + number_of_elements_per_line] = 0.25 * (
            three_part_ham[1][0] + three_part_ham[3][0]) ** 2

    # energy second translation on diff site
    array_with_energy_for_calculate[1][2 + number_of_elements_per_line] = 0.25 * (
            - 1 / math.sqrt(2) * two_part_ham[3][0] * (three_part_ham[1][0] + three_part_ham[3][0]) +
            two_part_ham[4][0] * (three_part_ham[0][0] + three_part_ham[2][0])) ** 2

    array_with_energy_for_calculate[2][2 + number_of_elements_per_line] = 0.25 * (
            - 1 / math.sqrt(2) * two_part_ham[3][0] * (three_part_ham[0][0] + three_part_ham[2][0]) +
            two_part_ham[5][0] * (three_part_ham[1][0] + three_part_ham[3][0])) ** 2

    # energy second translation on diff site
    array_with_energy_for_calculate[1][3 + number_of_elements_per_line] = 0.25 * (
            - 1 / math.sqrt(2) * two_part_ham[3][2] * (three_part_ham[1][0] + three_part_ham[3][0]) +
            two_part_ham[4][2] * (three_part_ham[0][0] + three_part_ham[2][0])) ** 2

    array_with_energy_for_calculate[2][3 + number_of_elements_per_line] = 0.25 * (
            - 1 / math.sqrt(2) * two_part_ham[3][2] * (three_part_ham[0][0] + three_part_ham[2][0]) +
            two_part_ham[5][2] * (three_part_ham[1][0] + three_part_ham[3][0])) ** 2

    # energy second translation on diff site
    array_with_energy_for_calculate[1][4 + number_of_elements_per_line] = 0.25 * (
            - 1 / math.sqrt(2) * two_part_ham[3][1] * (three_part_ham[1][0] + three_part_ham[3][0]) +
            two_part_ham[4][1] * (three_part_ham[0][0] + three_part_ham[2][0])) ** 2

    array_with_energy_for_calculate[2][4 + number_of_elements_per_line] = 0.25 * (
            - 1 / math.sqrt(2) * two_part_ham[3][1] * (three_part_ham[0][0] + three_part_ham[2][0]) +
            two_part_ham[5][1] * (three_part_ham[1][0] + three_part_ham[3][0])) ** 2

    number_of_elements_per_line += 5

    for i in range(2):
        array_with_energy_for_calculate[0][i + number_of_elements_per_line] = energ_of_state[
            i + number_of_elements_per_line - 18 * (syst_number - 1)]

        # energy first translation on diff site
    array_with_energy_for_calculate[1][0 + number_of_elements_per_line] = 0.5 * (
            three_part_ham[0][0] + three_part_ham[2][0]) ** 2
    array_with_energy_for_calculate[2][0 + number_of_elements_per_line] = 0.5 * (
            three_part_ham[1][0] + three_part_ham[3][0]) ** 2

    # energy first translation on diff site
    array_with_energy_for_calculate[1][1 + number_of_elements_per_line] = 0.5 * (
            three_part_ham[0][2] + three_part_ham[2][2]) ** 2
    array_with_energy_for_calculate[2][1 + number_of_elements_per_line] = 0.5 * (
            three_part_ham[1][2] + three_part_ham[3][2]) ** 2

    number_of_elements_per_line += 2

    return number_of_elements_per_line


# @njit(parallel=True)
def DOS_calc_for_inter_system(energ_array, zone_weight, number_of_sys, coulomb_potential):
    energy_distribution = energ_array[0][0] - 0.5
    while energy_distribution < energ_array[0][number_of_sys * 18 - 1] + 0.5:

        #element_index = binary_search_recursive(energ_array, energy_distribution, 0, len(energ_array[0]) - 1)
        dos = 0
        for array_element in range(0, len(energ_array[0]) - 1):
            dos = dos + (energ_array[1][array_element] + energ_array[2][array_element]) * \
                      approx_of_delta_function(0.0015, energ_array[0][array_element], energy_distribution)

        energy_distribution += 0.002

        output_to_file('DOS for system with interaction U = ' + str(coulomb_potential), str(dos / number_of_sys))


# Возможно не стоит подавать весь массив, если нужна только энергия и коэффициенты
# @njit(parallel=True)
def GIPR(energ_array, array_element, energy_distribution):
    gipr = ((energ_array[1][array_element] * approx_of_delta_function(0.0015, energ_array[0][array_element],
                                                                      energy_distribution)) ** 2 +
            (energ_array[2][array_element] * approx_of_delta_function(0.0015, energ_array[0][array_element],
                                                                      energy_distribution)) ** 2) / (
                   energ_array[1][array_element] * approx_of_delta_function(0.0015, energ_array[0][array_element],
                                                                            energy_distribution) +
                   energ_array[2][array_element] * approx_of_delta_function(0.0015, energ_array[0][array_element],
                                                                            energy_distribution)) ** 2
    return gipr


# допилить логику
# @njit(parallel=True)
def ensemble_averaged_GIPR(energ_array, zone_weight,number_of_sys,coulomb_potential):
    energy_distribution = energ_array[0][0] - 1
    while energy_distribution < energ_array[0][number_of_sys * 18 - 1] + 1:
        element_index = binary_search_recursive(energ_array, energy_distribution, 0, len(energ_array[0]) - 1)
        esemble_average_gipr = 0
        intermidate = 0

        for array_element in range(0, len(energ_array[0]) - 1):
                esemble_average_gipr = esemble_average_gipr + (
                        (energ_array[1][array_element] + energ_array[2][array_element]) * approx_of_delta_function(0.0015,
                                                                                           energ_array[0][array_element],
                                                                                           energy_distribution) *
                        GIPR(energ_array, array_element, energy_distribution) * approx_of_delta_function(0.0015,
                                                                                             energ_array[0][array_element],
                                                                                             energy_distribution))
                intermidate = intermidate + (energ_array[1][array_element] + energ_array[2][array_element]) * approx_of_delta_function(
                    0.0015, energ_array[0][array_element], energy_distribution) * (
                                  approx_of_delta_function(0.0015, energ_array[0][array_element], energy_distribution))

        esemble_average_gipr = esemble_average_gipr / intermidate

        output_to_file('IPR for system with interaction U = ' + str(coulomb_potential), str(esemble_average_gipr))
        energy_distribution += 0.002

def output_to_file(file_name, output_values):
    file = open(file_name + '.txt', 'a')

    file.write(output_values + '\n')

    file.close()


def energy_distribution(zone_weight, energ_array, number_of_sys,coulomb_potential):
    energy_distribution = energ_array[0][0] - 1
    while energy_distribution < energ_array[0][number_of_sys * 18 - 1] + 1:
        output_to_file('energy_distribution for system with interaction U = ' + str(coulomb_potential), str(energy_distribution))
        energy_distribution += 0.002


def quick_sort(mas, fst, lst):
    if fst >= lst:
        return mas

    i, j = fst, lst
    pivot = mas[0][int((fst + lst) / 2)]
    while i <= j:
        while mas[0][i] < pivot:
            i += 1
        while mas[0][j] > pivot:
            j -= 1
        if i <= j:
            mas[0][i], mas[0][j] = mas[0][j], mas[0][i]
            mas[1][i], mas[1][j] = mas[1][j], mas[1][i]
            mas[2][i], mas[2][j] = mas[2][j], mas[2][i]
            i, j = i + 1, j - 1

    quick_sort(mas, fst, j)
    quick_sort(mas, i, lst)

def array_out(array):

    for i in range(len(array[0])):
        output_to_file("array_with_energy", str(array[0][i]))

def binary_search_recursive(array, element, start, end):
    if start > end:
        return -1

    mid = (start + end) // 2
    if int(array[0][mid] * 1000)/1000 - 0.05 <= element <= int(array[0][mid] * 1000)/1000 + 0.05:
        return mid

    if element < int(array[0][mid]*1000)/1000 - 0.05:
        return binary_search_recursive(array, element, start, mid-1)
    if element > int(array[0][mid] * 1000) / 1000 + 0.05:
        return binary_search_recursive(array, element, mid+1, end)
