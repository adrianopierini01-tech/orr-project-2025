from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from ase.io import read
from fairchem.core import pretrained_mlip, FAIRChemCalculator
import numpy as np

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calculator = FAIRChemCalculator(predictor, task_name="oc20")

# Create the metal slab
slab = fcc111("Pt", (3, 3, 4), vacuum=12, periodic=True)
slab.center(vacuum=12, axis=2)

# Freeze the bottom Pt layers
fix = FixAtoms(indices=range(0,18))
slab.set_constraint(fix)

# Add the adsorbate
adsorbate = "O"
add_adsorbate(slab, adsorbate, 2.0, "fcc")

# Create a water monolayer
water = read("molecules/h2o.xyz")
water_layer = []
top_sites = [28,31,32,33,35] # surface site where water molecules will be placed
for index in top_sites:
    w = water.copy()
    new_position = slab[index].position + [0.0, 0.0 ,2.0] # translate H2O on top of sites
    w.translate(new_position)
    w.euler_rotate( # random rotation of H2O
                    ... # HINT # (phi, theta, psi, center) use np.random.random
                   center=w.get_center_of_mass(),
                   )
    water_layer.append(w)

# Add the water molecules to the slab+adsorbate
water_atoms = Atoms()
for w in water_layer:
    water_atoms = water_atoms + w
system = slab + water_atoms

# Write the initial geometry to file
system.write("init.pt-o.xyz", format="xyz")

# Associate the calculator
system.calc = calculator

# Set up LBFGS dynamics object
opt = LBFGS(atoms=system)
opt.run(0.01, 250)

# Write the final geometry to file
system.write("final.pt-o.xyz", format="xyz")
