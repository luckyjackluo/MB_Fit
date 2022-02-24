import warnings
# absolute module imports
from mbfit.utils import constants, SettingsReader, files
from mbfit.exceptions import NoEnergiesError, NoOptimizedEnergyError, MultipleOptimizedEnergiesError, NoEnergyInRangeError

# local module imports
from .database import Database
from mbfit.utils import system


def generate_1b_training_set(settings_path, database_config_path, training_set_path, molecule_name, method, basis, cp, *tags, e_min=0, e_max=float('inf')):
    """
    Writes a 1b training set to the given file from the calculated energies in a database.

    ***deprecated, please use generate_training_set instead***

    Args:
        settings_path       - Local path to the ".ini" file with all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path   - Local path to file to write training set to.
        molecule_name       - The name of the molecule to generate a training set for.
        method              - Use energies calculated with this method. Use % for any method.
        basis               - Use energies calculated with this basis. Use % for any basis.
        cp                  - Use energies calculated with this cp. Use 0 for False, 1 for True, or % for any cp.
        tags                - Use energies marked with one or more of these tags. Use % for any tag.
        e_min               - The minimum (inclusive) energy of any configuration to include in the training set.
        e_max               - The maximum (exclusive) energy of any configuration to include in the training set.

    Return:
        None.
    """

    warnings.warn("generate_1b_training_set() has been deprecated. Please use generate_training_set() instead.",
            DeprecationWarning)

    settings = SettingsReader(settings_path)

    names = settings.get("molecule", "names").split(",")
    SMILES = settings.get("molecule", "SMILES").split(",")
    
    # open the database
    with Database(database_config_path) as database:

        print("Creating a fitting input file from database into file {}".format(training_set_path))

        #intializing a counter
        count_configs = 0

        with open(files.init_file(training_set_path, files.OverwriteMethod.get_from_settings(settings)), "w") as output:

            for molecule, deformation_energy_hartrees in database.get_1B_training_set(molecule_name, names, SMILES, method, basis, cp, *tags):

                deformation_energy_kcalmol = deformation_energy_hartrees * constants.au_to_kcal # Converts Hartree to kcal/mol
                if deformation_energy_kcalmol < e_max and deformation_energy_kcalmol - e_min > -0.000000000001:
                    # write the number of atoms to the output file
                    output.write(str(molecule.get_num_atoms()) + "\n")
                    output.write(str(deformation_energy_kcalmol) + " ")
                    output.write("\n")

                    # write the molecule's atoms' coordinates to the xyz file
                    output.write(molecule.to_xyz() + "\n")

                    # increment the counter
                    count_configs += 1

            if count_configs == 0:
                raise Exception

            print("Generated training set with " + str(count_configs) + " Configurations.")


def generate_2b_training_set(settings_path, database_config_path, training_set_path, molecule_name, method, basis,
        cp, *tags, e_bind_max=float('inf'), e_mon_max=float('inf')):
    """"
    Creates a 2b training set file from the calculated energies in a database.

    ***deprecated, please use generate_training_set instead***

    Args:
        settings_path       - Local path to the ".ini" file with all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path   - Local path to file to write training set to.
        molecule_name       - The name of this dimer.
        method              - Use energies calculated with this method. Use % for any method.
        basis               - Use energies calculated with this basis. Use % for any basis.
        cp                  - Use energies calculated with this cp. Use 0 for False, 1 for True, or % for any cp.
        tags                - Use energies marked with at least one of these tags. Use % for any tag.
        e_bind_max          - Maximum binding energy allowed
        e_mon_max           - Maximum monomer deformation energy allowed

    Return:
        None.
    """

    warnings.warn("generate_2b_training_set() has been depricated. Please use generate_training_set() instead.",
            DeprecationWarning)

    settings = SettingsReader(settings_path)

    names = settings.get("molecule", "names").split(",")
    SMILES = settings.get("molecule", "SMILES").split(",")
    
    # open the database
    with Database(database_config_path) as database:

        print("Creating a fitting input file from database into file {}".format(training_set_path))

        # initializing a counter
        count_configs = 0

        with open(files.init_file(training_set_path, files.OverwriteMethod.get_from_settings(settings)), "w") as output:
            for molecule, binding_energy, interaction_energy, monomer1_energy_deformation, monomer2_energy_deformation in database.get_2B_training_set(molecule_name, names, SMILES, method, basis, cp, *tags):

                binding_energy *= constants.au_to_kcal
                interaction_energy *= constants.au_to_kcal
                monomer1_energy_deformation *= constants.au_to_kcal
                monomer2_energy_deformation *= constants.au_to_kcal

                # write the configuration only if is below the thresholds
                if (binding_energy < e_bind_max) and (monomer1_energy_deformation < e_mon_max) and (
                        monomer2_energy_deformation < e_mon_max):

                    # write the number of atoms to the output file
                    output.write(str(molecule.get_num_atoms()) + "\n")

                    output.write("{} {} {} {}".format(binding_energy, interaction_energy, monomer1_energy_deformation,
                                                      monomer2_energy_deformation))

                    output.write("\n")

                    # write the molecule's atoms' coordinates to the xyz file
                    output.write(molecule.to_xyz() + "\n")

                    # increment the counter
                    count_configs += 1

        if count_configs == 0:
            raise Exception

        print("Generated training set with " + str(count_configs) + " Configurations.")


def generate_training_set(settings_path, database_config_path, training_set_path, method, basis,
        cp, *tags, e_bind_min=-float('inf'), e_bind_max=float('inf'), e_mon_min=-float('inf'), e_mon_max=float('inf'),
        deprecated_fitcode=False):
    """"
    Creates a training set file from the calculated energies in a database.

    Args:
        settings_path       - Local path to the ".ini" file with all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path   - Local path to file to write training set to.
        method              - Use energies calculated with this method. Use % for any method.
        basis               - Use energies calculated with this basis. Use % for any basis.
        cp                  - Use energies calculated with this cp. Use 0 for False, 1 for True, or % for any cp.
        tags                - Use energies marked with at least one of these tags. Use % for any tag.
        e_bind_min          - Minimum binding energy allowed, inclusive.
        e_bind_max          - Maximum binding energy allowed, exclusive.
        e_mon_max           - Minimum monomer deformation energy allowed, inclusive.
        e_mon_max           - Maximum monomer deformation energy allowed, exclusive.
        deprecated_fitcode  - Is this function being called to be used with the deprecated fitcode?
                The output of the 1b and 2b training sets will be different.

    Return:
        None.
    """

    settings = SettingsReader(settings_path)

    names = settings.get("molecule", "names").split(",")
    SMILES = settings.get("molecule", "SMILES").split(",")

    # open the database
    with Database(database_config_path) as database:

        training_set_size = database.get_training_set_size(names, method, basis, cp, *tags)
        system.format_print("Creating a fitting input file from database into file {} with up to {} geometries.".format(training_set_path, training_set_size),
                bold=True, color=system.Color.YELLOW)

        # initializing a counter
        count_configs = 0
        filtered_configs = 0

        with open(files.init_file(training_set_path, files.OverwriteMethod.get_from_settings(settings)), "w") as output:
            for molecule, binding_energy, nb_energy, deformation_energies in database.get_training_set(names, SMILES, method, basis, cp, *tags):

                binding_energy *= constants.au_to_kcal
                nb_energy *= constants.au_to_kcal
                deformation_energies = [d * constants.au_to_kcal for d in deformation_energies]

                # skip this config if the binding energy is >= the maximum.
                if binding_energy < e_bind_min or binding_energy >= e_bind_max:
                    filtered_configs += 1
                    continue

                # skip this config if any binding energy is >= the maximum.
                if any([d < e_mon_min or d >= e_mon_max for d in deformation_energies]):
                    filtered_configs += 1
                    continue

                # write the number of atoms to the output file
                output.write(str(molecule.get_num_atoms()) + "\n")

                if deprecated_fitcode:
                    if molecule.get_num_fragments() == 1:
                        output.write("{}".format(binding_energy))
                    elif molecule.get_num_fragments() == 2:
                        output.write("{} {} {} {}".format(binding_energy, nb_energy, deformation_energies[0], deformation_energies[1]))
                    else:
                        output.write("{} {}".format(binding_energy, nb_energy))
                else:
                    output.write("{} {}".format(binding_energy, nb_energy))

                output.write("\n")

                # write the molecule's atoms' coordinates to the xyz file
                output.write(molecule.to_xyz() + "\n")

                # increment the counter
                count_configs += 1

                if count_configs + filtered_configs % 100 == 0:
                    system.format_print("Considered {} geometries so far. Included {} and filtered {} out.".format(count_configs + filtered_configs, count_configs, filtered_configs),
                            italics=True)

        if count_configs == 0:
            raise Exception

        system.format_print("Generated training set with {} configurations. {} configurations filtered out due to binding or deformation energies outside of specified range.".format(count_configs, filtered_configs),
                bold=True, color=system.Color.GREEN)
