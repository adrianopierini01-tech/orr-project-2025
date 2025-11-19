import numpy as np

# -----------------------------
# USER SETTINGS
# -----------------------------

# Input files
file_OOH = "../4_eads/eads.ooh.txt"
file_O = "../4_eads/eads.o.txt"

# Output file
file_output = "delta_G1.txt"

# Potential bias
U = ... # HINT # potential in Volts

# Hard-coded reference energies
G_H2O = 0.00
correct_OOH = -0.22
correct_O = +0.00

# -----------------------------
# Delta_G for this reaction step
# -----------------------------

def deltaG(E_OOH, E_O):
    return (E_O + correct_O) + G_H2O - (E_OOH + correct_OOH) + U

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
data_OOH = read_energy_file(file_OOH)

# Ensure both files contain the same labels
labels = sorted(set(data_O.keys()))

# -----------------------------
# WRITE output files
# -----------------------------

with open(file_output, "w") as f:

    f.write("Metal\tdG_reaction_eV\n")

    for label in labels:
        E_O = data_O[label]
        E_OOH = data_OOH[label]

        dG = deltaG(E_OOH, E_O)

        f.write(f"{label}\t{dG:.6f}\n")

print(f"Done. Reaction free energies written to: {file_output}")
