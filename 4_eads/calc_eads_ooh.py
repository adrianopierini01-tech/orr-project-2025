import numpy as np

# -----------------------------
# USER SETTINGS
# -----------------------------

# Input files
file_OOH = "../3_fcc_metals/etot.metal-ooh.txt"
file_slab = "../3_fcc_metals/etot.metal.txt"

# Output file
file_output = "eads.ooh.txt"

# Hard-coded reference energies
E_H2 = ... # HINT # gas-phase energy
E_H2O = ... # HINT # gas-phase energy

# -----------------------------
# GENERAL FORMULAS
# Replace these with any formula you need
# -----------------------------

def delta_Eads(E_slab, E_OOH):
    return ( E_OOH + 1.5*E_H2 ) - ( E_slab + 2.0*E_H2O )

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
data_slab = read_energy_file(file_slab)

# Ensure both files contain the same labels
labels = sorted(set(data_OOH.keys()) & set(data_slab.keys()))

# -----------------------------
# WRITE output files
# -----------------------------

with open(file_output, "w") as f:

    f.write("Metal\tEnergy_ads_eV\n")

    for label in labels:
        E_OOH = data_OOH[label]
        E_slab= data_slab[label]

        Eads = delta_Eads(E_slab, E_OOH)

        f.write(f"{label}\t{Eads:.6f}\n")

print(f"Done. Reaction free energies written to: {file_output}")
