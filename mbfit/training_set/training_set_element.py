class TrainingSetElement:

    def __init__(self, molecule, **energies):
        self.molecule = molecule

        self.energies = energies

    def get_molecule(self):
        return self.molecule

    def has_energy(self, energy_name):
        return energy_name in self.energies.keys()

    def get_energy(self, energy_name):
        return self.energies.get(energy_name)

    def add_energy(self, energy_name, energy):
        self.energies[energy_name] = energy

    def __str__(self):
        return self.molecule.to_xyz() + "\n" + str(self.energies)
