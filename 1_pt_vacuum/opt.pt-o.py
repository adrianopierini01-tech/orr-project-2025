from ase.build import fcc111, add_adsorbate
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from fairchem.core import pretrained_mlip, FAIRChemCalculator

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calculator = FAIRChemCalculator(predictor, task_name="oc20")

# Create the metal slab
slab = fcc111("Pt", (3, 3, 4), vacuum=12, periodic=True)
slab.center(vacuum=12, axis=2)

# Freeze the bottom Pt layers
fix = FixAtoms(indices=range(0,18))
slab.set_constraint(fix)

# add the adsorbate
adsorbate = "O"
add_adsorbate(...)   # HINT # ( surface, adsorbate, dist, "adsorp_site" )

# Associate the calculator
slab.calc = calculator

# Set up LBFGS dynamics object
opt = LBFGS(atoms=slab)
opt.run(0.01, 250)

# Write optimized structure
slab.write("final.pt-o.xyz", format="xyz")
