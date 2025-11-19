import numpy as np

# -----------------------------
# USER SETTINGS
# -----------------------------

# Input files
file_OH = "../4_eads/eads.oh.txt"
file_O = "../4_eads/eads.o.txt"

# Output file
file_output = "delta_G2.txt"

# Potential bias
U = ... # HINT # potential in Volts

# Hard-coded reference energies
correct_OH = -0.30
correct_O = +0.00

# -----------------------------
# Delta_G for this reaction step
# -----------------------------

def deltaG(E_O, E_OH):
    return (E_OH + correct_OH) - (E_O + correct_O) + U

# -----------------------------
# READ input files
# Expect format: label   energy
# -----------------------------

def read_energy_file(filename):
    data = {}
    with open(filename, "r") as f:
        for line in f.readlines()[1:]:
            if not line.strip():
                continue
            label, val = line.split()
            data[label] = float(val)
    return data

data_O = read_energy_file(file_O)
data_OH = read_energy_file(file_OH)

# Ensure both files contain the same labels
labels = sorted(set(data_O.keys()))

# -----------------------------
# WRITE output files
# -----------------------------

with open(file_output, "w") as f:

    f.write("Metal\tdG_reaction_eV\n")

    for label in labels:
        E_O = data_O[label]
        E_OH = data_OH[label]

        dG = deltaG(E_O, E_OH)

        f.write(f"{label}\t{dG:.6f}\n")

print(f"Done. Reaction free energies written to: {file_output}")
