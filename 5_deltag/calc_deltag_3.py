import numpy as np

# -----------------------------
# USER SETTINGS
# -----------------------------

# Input files
file_OH = "../4_eads/eads.oh.txt"

# Output file
file_output = "delta_G3.txt"

# Potential bias
U = ... # HINT # potential in Volts

# Hard-coded reference energies
G_H2O = 0.00
correct_OH = -0.30

# -------------------------------
# Delta_G for this reaction step
# -------------------------------

def deltaG(E_OH):
    return G_H2O - (E_OH + correct_OH) + U

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

data_OH = read_energy_file(file_OH)

# Ensure both files contain the same labels
labels = sorted(set(data_OH.keys()))

# -----------------------------
# WRITE output files
# -----------------------------

with open(file_output, "w") as f:

    f.write("Metal\tdG_reaction_eV\n")

    for label in labels:
        E_OH = data_OH[label]

        dG = deltaG(E_OH)

        f.write(f"{label}\t{dG:.6f}\n")

print(f"Done. Reaction free energies written to: {file_output}")
