from ase.build import fcc111, add_adsorbate
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from ase.io import read
from fairchem.core import pretrained_mlip, FAIRChemCalculator

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calculator = FAIRChemCalculator(predictor, task_name="oc20")

# Create the metal slab
...
slab.center(vacuum=12, axis=2)

# Freeze the bottom Pt layers
fix = FixAtoms(indices=range(0,18))
slab.set_constraint(fix)

# add the adsorbate
adsorbate = read("molecules/ooh.xyz")
add_adsorbate(slab, adsorbate, 2.0, "ontop")

# Associate the calculator
slab.calc = calculator

# Set up LBFGS dynamics object
opt = LBFGS(atoms=slab)
opt.run(0.01, 250)

# Write optimized structure
slab.write("final.alloy-ooh.xyz", format="xyz")
