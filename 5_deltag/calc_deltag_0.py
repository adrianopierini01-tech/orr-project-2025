import numpy as np

# -----------------------------
# USER SETTINGS
# -----------------------------

# Input files
file_OOH = "../4_eads/eads.ooh.txt"

# Output file
file_output = "delta_G0.txt"

# Potential bias
U = ... # HINT # potential in Volts

# Hard-coded reference energies
G_O2 = 4.92
correct_OOH = -0.22

# -----------------------------
# Delta_G for this reaction step
# -----------------------------

def deltaG(E_OOH):
    return ( E_OOH + correct_OOH ) - G_O2 + U

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

data_OOH = read_energy_file(file_OOH)

# Ensure both files contain the same labels
labels = sorted(set(data_OOH.keys()))

# -----------------------------
# WRITE output files
# -----------------------------

with open(file_output, "w") as f:

    f.write("Metal\tdG_reaction_eV\n")

    for label in labels:
        E_OOH = data_OOH[label]

        dG = deltaG(E_OOH)

        f.write(f"{label}\t{dG:.6f}\n")

print(f"Done. Reaction free energies written to: {file_output}")
