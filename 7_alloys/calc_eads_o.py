import numpy as np

# -----------------------------
# USER SETTINGS
# -----------------------------

# Input files
file_O = "./etot.alloy-o.txt"
file_slab = "./etot.alloy.txt"

# Output file
file_output = "eads.o.txt"

# Hard-coded reference energies
E_H2 = ... # HINT # gas-phase energy
E_H2O = ... # HINT # gas-phase energy

# -----------------------------
# Calculate E_ads
# -----------------------------

def delta_Eads(E_slab, E_O):
    return ( E_O + E_H2 ) - ( E_slab + E_H2O )

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

data_O = read_energy_file(file_I1)
data_slab = read_energy_file(file_slab)

# Ensure both files contain the same labels
labels = sorted(set(data_O.keys()) & set(data_slab.keys()))

# -----------------------------
# WRITE output files
# -----------------------------

with open(file_output, "w") as f:

    f.write("Metal\tEnergy_ads_eV\n")

    for label in labels:
        E_O = data_I1[label]
        E_slab= data_slab[label]

        Eads = delta_Eads(E_slab, E_O)

        f.write(f"{label}\t{Eads:.6f}\n")

print(f"Done. Reaction free energies written to: {file_output}")
