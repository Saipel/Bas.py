# @njit(parallel=True)
def DOS_calc_for_inter_system(energ_array, zone_weight, number_of_sys, coulomb_potential):
    energy_distribution = energ_array[0][0] - 0.5
    while energy_distribution < energ_array[0][number_of_sys * 18 - 1] + 0.5:

        element_index = binary_search_recursive(energ_array, energy_distribution, 0, len(energ_array[0]) - 1)
        dos = 0
        if element_index < 250_001:
            for array_element in range(0, element_index + 250_000):
                dos = dos + (energ_array[1][array_element] + energ_array[2][array_element]) * \
                      approx_of_delta_function(0.0015, energ_array[0][array_element], energy_distribution)

        if 250_001 <= element_index < len(energ_array[0]) - 250_000:
            for array_element in range(element_index - 250_000, element_index + 250_000):
                dos = dos + (energ_array[1][array_element] + energ_array[2][array_element]) * \
                      approx_of_delta_function(0.0015, energ_array[0][array_element], energy_distribution)

        if element_index >= len(energ_array[0]) - 250_000:
            for array_element in range(element_index - 250_000, len(energ_array[0]) - 1):
                dos = dos + (energ_array[1][array_element] + energ_array[2][array_element]) * \
                      approx_of_delta_function(0.0015, energ_array[0][array_element], energy_distribution)

        energy_distribution += 0.002

        output_to_file('DOS for system with interaction U = ' + str(coulomb_potential), str(dos / number_of_sys))


# допилить логику
# @njit(parallel=True)
def ensemble_averaged_GIPR(energ_array, zone_weight,number_of_sys,coulomb_potential):
    energy_distribution = energ_array[0][0] - 1
    while energy_distribution < energ_array[0][number_of_sys * 18 - 1] + 1:
        element_index = binary_search_recursive(energ_array, energy_distribution, 0, len(energ_array[0]) - 1)
        esemble_average_gipr = 0
        intermidate = 0

        if element_index < 250_000:
            for array_element in range(0, element_index + 250_000):
                esemble_average_gipr = esemble_average_gipr + (
                        (energ_array[1][array_element] + energ_array[2][array_element]) * approx_of_delta_function(0.0015, energ_array[0][array_element],
                                                                                           energy_distribution) *
                        GIPR(energ_array, array_element, energy_distribution) * approx_of_delta_function(0.0015, energ_array[0][array_element],
                                                                                             energy_distribution))
                intermidate = intermidate + (energ_array[1][array_element] + energ_array[2][array_element]) * approx_of_delta_function(0.0015,
                                                                                                               energ_array[
                                                                                                                   0][
                                                                                                                   array_element],
                                                                                                               energy_distribution) * (
                                  approx_of_delta_function(0.0015, energ_array[0][array_element], energy_distribution))
            esemble_average_gipr = esemble_average_gipr / intermidate

        if 250_000 <= element_index < len(energ_array[0]) - 250_000:
            for array_element in range(element_index - 250_000, element_index + 250_000):
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

        if element_index >= len(energ_array[0]) - 250_000:
            for array_element in range(element_index - 250_000, len(energ_array[0]) - 1):
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