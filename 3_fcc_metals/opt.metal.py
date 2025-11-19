from ase.build import fcc111
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from fairchem.core import pretrained_mlip, FAIRChemCalculator

# -------- Settings --------
metals = ["Ag","Au","Cu","Ir","Pt","Pd","Rh"]
fmax = 0.01
steps = 250
outfile = "etot.metal.txt"

# -------- ML potential --------
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calculator = FAIRChemCalculator(predictor, task_name="oc20")

# ----- Optimization routine ----
def optimize_slab(metal):

    print(f"\n=== Optimizing {metal}(111) ===")

    # Build slab
    slab = fcc111(metal, (3, 3, 4), vacuum=12, periodic=True)
    slab.center(vacuum=12, axis=2)

    # Freeze bottom layers
    fix = FixAtoms(indices=range(0,18))
    slab.set_constraint(fix)

    # Attach calculator
    slab.calc = calculator

    # Optimize
    opt = LBFGS(slab)
    opt.run(fmax=fmax, steps=steps)

    return slab.get_potential_energy()


# -------- Run optimizations --------
with open(outfile, "w") as f:
    f.write("Metal\tFinal_Energy_eV\n")

    ... # HINT # iterate over metals
        ... # HINT # call the optimization routine

        print(f"Final energy for {metal}: {energy:.6f} eV")
        f.write(f"{metal}\t{energy:.6f}\n")

print(f"\nDone. Energies saved to {outfile}")
