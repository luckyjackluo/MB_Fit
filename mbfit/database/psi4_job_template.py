import psi4
from psi4.driver.qcdb.exceptions import QcdbException
import subprocess, os

"""
This is not a runnable python file.

It is a template that database_job_maker.py uses to make psi4 jobs
"""

def execute_job():
    whole_molecule = "{whole_molecule}"
    charges = "{charges}"
    spins = "{spins}"
    symmetries = "{symmetries}"
    SMILES = "{SMILES}"
    names = "{names}"
    atom_counts = "{atom_counts}"
    total_atoms = "{total_atoms}"
    molecule = "{molecule}"
    frag_indices = "{frag_indices}"
    method = "{method}"
    basis = "{basis}"
    cp = "{cp}"
    use_cp = "{use_cp}"
    number_of_threads = {num_threads}
    memory = "{memory}"
    total_charge = "{total_charge}"
    total_spin = "{total_spin}"
    job_hash = "{job_hash}"
    qm_options = {qm_options}

    try:
        max_threads = int(subprocess.check_output(["grep", "-c", "cores", "/proc/cpuinfo"]))
        print("Maximum threads: {format}".format(max_threads))
    except:
        print("Error detecting number of cores. \n Maxium threads: 1")
        max_threads = 1

    if number_of_threads > max_threads:
        print("Input number of threads ({format}) greater than max threads ({format}), limiting number of threads to max threads".format(number_of_threads, max_threads))
        number_of_threads = max_threads

    print("Running Job")
    print("Molcule: {format}".format(molecule))
    print("Method: {format}".format(method))
    print("Basis: {format}".format(basis))
    print("Threads {format}".format(number_of_threads))
    print("Memory: {format}".format(memory))

    i = 8
    job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "job_{format}".format(job_hash[:i]))

    while os.path.exists(job_dir):
        i += 1

        job_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "job_{format}".format(job_hash[:i]))

    if os.path.isdir(job_dir + "_done"):
        print("Job " + job_dir.split('/')[-1] + " is already done. Skipping.")
        return

    os.mkdir(job_dir)
    output_path = job_dir + "/output.ini"
    log_path = job_dir + "/output.log"

    psi4.set_options(qm_options)

    # psi4 specific stuff!

    psi4.core.set_output_file(log_path, False)
    psi4.set_memory(memory)
    psi4_input_geo = "\n" + str(total_charge) + " " + str(total_spin) + "\n" + molecule
    psi4.geometry(psi4_input_geo)
    psi4.set_num_threads(number_of_threads)

    try:
        energy = psi4.energy("{format}/{format}".format(method, basis))
        print("Energy: {format}".format(energy))
        success = True
    except (ValueError, SystemError):
        success = False
        print("Something went wrong...")
    except QcdbException:
        success = False
        print("The calculation failed.")
    except RuntimeError:
        success = False
        print("The calculation failed. Iterations *probably* did not converge.")

    # end psi4 specific stuff

    with open(output_path, "w") as out_file:
        out_file.write("[molecule]\n")
        out_file.write("xyz = {format}\n\n{format}".format(total_atoms, whole_molecule).replace("\n", "\n ") + "\n")
        out_file.write("atom_counts = {format}\n".format(atom_counts))
        out_file.write("charges = {format}\n".format(charges))
        out_file.write("spins = {format}\n".format(spins))
        out_file.write("symmetries = {format}\n".format(symmetries).replace("'", ""))
        out_file.write("SMILES = {format}\n".format(SMILES).replace("'", ""))
        out_file.write("names = {format}\n".format(names).replace("'",""))
        out_file.write("method = {format}\n".format(method))
        out_file.write("basis = {format}\n".format(basis))
        out_file.write("cp = {format}\n".format(cp))
        out_file.write("use_cp = {format}\n".format(use_cp))
        out_file.write("frag_indices = {format}\n".format(frag_indices))
        out_file.write("job_hash = {format}\n".format(job_hash))

        if success:
            out_file.write("energy = {format}\n".format(energy))


if __name__ == "__main__":
    execute_job()
