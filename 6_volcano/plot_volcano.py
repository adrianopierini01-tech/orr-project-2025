# ==============================================================
# USER INPUT
# ==============================================================

import numpy as np
import matplotlib.pyplot as plt

# Define functions for Ulim for each step
def Ulim_0(E_ads_OOH):
    return 5.14 - E_ads_OOH

def Ulim_1(E_ads_OH, E_ads_OOH): # A scaling relation dG(O) = 2*dG(OH) was used
    return E_ads_OOH - 2*E_ads_OH + 0.38

def Ulim_2(E_ads_OH): # A scaling relation dG(O) = 2*dG(OH) was used
    return E_ads_OH - 0.30

def Ulim_3(E_ads_OH):
    return E_ads_OH - 0.30

# Define U_lim as a function of E_ads_O and E_ads_OH
def Ulim_all(E_ads_OOH, E_ads_OH):
    u0 = Ulim_0(E_ads_OOH)
    u1 = Ulim_1(E_ads_OH, E_ads_OOH)
    u2 = Ulim_2(E_ads_OH)
    u3 = Ulim_3(E_ads_OH)
    return min(u0, u1, u2, u3)

# Create a grid of E_ads_O and E_ads_OH values
E_ads_OH_values = np.linspace(0, 2.5, 100)
E_ads_OOH_values = np.linspace(2.5, 5, 100)
E_ads_OH_grid, E_ads_OOH_grid = np.meshgrid(E_ads_OH_values, E_ads_OOH_values)

# Calculate U_lim for each point in the grid
Ulim_values = np.zeros_like(E_ads_OH_grid)
for i in range(E_ads_OH_grid.shape[0]):
    for j in range(E_ads_OH_grid.shape[1]):
        Ulim_values[i, j] = Ulim_all(E_ads_OOH_grid[i, j], E_ads_OH_grid[i, j])

# Representative points for real systems (example values)
real_systems = [
    {"name": "Pt", "E_ads_O": 1.58, "E_ads_OOH": 4.05, "E_ads_OH": 1.10},
    {"name": "Pd", "E_ads_O": 1.38, "E_ads_OOH": 4.06, "E_ads_OH": 0.90},
    {"name": "Au", "E_ads_O": 2.77, "E_ads_OOH": 4.65, "E_ads_OH": 1.53},
    {"name": "Ag", "E_ads_O": 2.22, "E_ads_OOH": 4.17, "E_ads_OH": 0.84},
    {"name": "Rh", "E_ads_O": 0.54, "E_ads_OOH": 3.77, "E_ads_OH": 0.42},
    {"name": "Ir", "E_ads_O": 0.99, "E_ads_OOH": 3.85, "E_ads_OH": 0.78},
    {"name": "Cu", "E_ads_O": 1.05, "E_ads_OOH": 3.71, "E_ads_OH": 0.36}
]

# Plot the volcano plot as a color map
plt.figure(figsize=(8, 6))
contour = plt.contourf(E_ads_OH_grid, E_ads_OOH_grid, Ulim_values, levels=20, cmap='viridis')
plt.colorbar(contour, label='U_lim (V)')
plt.xlabel('E_ads_OH (eV)')
plt.ylabel('E_ads_OOH (eV)')
plt.title('Volcano Plot: U_lim as a Function of E_ads_OH and E_ads_O')

# Plot representative points
for system in real_systems:
    plt.scatter(system["E_ads_OH"], system["E_ads_OOH"], color='red', label=system["name"])
    plt.text(system["E_ads_OH"], system["E_ads_OOH"], system["name"], color='black', fontsize=12, ha='center')

# Save the plot
plt.savefig('volcano.png')

# Show the plot
plt.show()
