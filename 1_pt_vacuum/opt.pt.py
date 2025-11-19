from ase.build import fcc111
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from ase.io import read
from fairchem.core import pretrained_mlip, FAIRChemCalculator

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calculator = FAIRChemCalculator(predictor, task_name="oc20")

# Create the metal slab
slab = fcc111( . . .)  # HINT # ("element", (Nx,Ny,Nz), d_vacuum, periodic=True/False)
slab.center(vacuum=12, axis=2)

# Freeze the bottom Pt layers
fix = FixAtoms(indices=range(0,18))
slab.set_constraint(fix)

# Associate the calculator
slab.calc = calculator

# Set up LBFGS dynamics object
opt = LBFGS(atoms=slab)
opt.run(0.01, 250)

# Write optimized structure
slab.write("final.pt.xyz", format="xyz")
