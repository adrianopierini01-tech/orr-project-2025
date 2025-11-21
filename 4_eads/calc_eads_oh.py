import numpy as np

# -----------------------------
# USER SETTINGS
# -----------------------------

# Input files
file_OH = "../3_fcc_metals/etot.metal-oh.txt"
file_slab = "../3_fcc_metals/etot.metal.txt"

# Output file
file_output = "eads.oh.txt"

# Hard-coded reference energies
E_H2 = ... # HINT # gas-phase energy
E_H2O = ... # HINT # gas-phase energy

# -----------------------------
# Calculate E_ads
# -----------------------------

def delta_Eads(E_slab, E_OH):
    return ( E_OH + 0.5*E_H2 ) - ( E_slab + E_H2O )

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
data_slab = read_energy_file(file_slab)

# Ensure both files contain the same labels
labels = sorted(set(data_OH.keys()) & set(data_slab.keys()))

# -----------------------------
# WRITE output files
# -----------------------------

with open(file_output, "w") as f:

    f.write("Metal\tEnergy_ads_eV\n")

    for label in labels:
        E_OH = data_OH[label]
        E_slab= data_slab[label]

        Eads = delta_Eads(E_slab, E_OH)

        f.write(f"{label}\t{Eads:.6f}\n")

print(f"Done. Reaction free energies written to: {file_output}")
