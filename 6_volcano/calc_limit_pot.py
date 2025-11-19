# ==============================================================
# USER INPUT
# ==============================================================

def read_energy_file(filename):
    data = {}
    with open(filename, "r") as f:
        for line in f.readlines()[1:]:
            if not line.strip():
                continue
            label, val = line.split()
            data[label] = float(val)
    return data

# Define functions for Ulim for each step
def Ulim_0(E_ads_OOH):
    return 5.14 - E_ads_OOH

def Ulim_1(E_ads_O, E_ads_OOH):
    return E_ads_OOH - E_ads_O - 0.22

def Ulim_2(E_ads_OH, E_ads_O):
    return E_ads_O - E_ads_OH + 0.30

def Ulim_3(E_ads_OH):
    return E_ads_OH - 0.30

# Define U_lim as a function of E_ads_O and E_ads_OH
def Ulim_all(E_ads_OOH, E_ads_O, E_ads_OH):
    u0 = Ulim_0(E_ads_OOH)
    u1 = Ulim_1(E_ads_O, E_ads_OOH)
    u2 = Ulim_2(E_ads_OH, E_ads_O)
    u3 = Ulim_3(E_ads_OH)
    return min(u0, u1, u2, u3)


# ==============================================================
# MAIN EXECUTION
# ==============================================================

data_OOH = read_energy_file("../4_eads/eads.ooh.txt")
data_O = read_energy_file("../4_eads/eads.o.txt")
data_OH = read_energy_file("../4_eads/eads.oh.txt")

file_output = "limiting_potentials.txt"

labels = sorted(set(data_OOH.keys()))

# -----------------------------
# WRITE output files
# -----------------------------

with open(file_output, "w") as f:

    f.write("Metal\tLim_Potential_V\n")

    for label in labels:
        E_ads_OOH = data_OOH[label]
        E_ads_O = data_O[label]
        E_ads_OH = data_OH[label]

        U = Ulim_all(E_ads_OOH, E_ads_O, E_ads_OH)

        f.write(f"{label}\t{U:.4f}\n")

print(f"Done. Reaction free energies written to: {file_output}")
