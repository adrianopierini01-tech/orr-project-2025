from ase.build import fcc111, add_adsorbate
from ase.optimize import LBFGS
from ase.constraints import FixAtoms
from fairchem.core import pretrained_mlip, FAIRChemCalculator

# -------- Settings --------
metals = ["Ag","Au","Cu","Ir","Ni","Pt","Pd","Rh"]
adsorbate = "O"
ads_height = 2.0
ads_site = "fcc"
fmax = 0.01
steps = 250
outfile = "etot.metal-o.txt"

# -------- ML potential --------
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calculator = FAIRChemCalculator(predictor, task_name="oc20")

# -------- Run optimizations --------
with open(outfile, "w") as f:
    f.write("Metal\tFinal_Energy_eV\n")

    for metal in metals:
        print(f"\n=== Optimizing {adsorbate} on {metal}(111) ===")

        # Build slab
        slab = fcc111(metal, (3, 3, 4), vacuum=12, periodic=True)
        slab.center(vacuum=12, axis=2)

        # Freeze bottom layers
        fix = FixAtoms(indices=range(0,18))
        slab.set_constraint(fix)

        # Add adsorbate
        add_adsorbate(slab, adsorbate, ads_height, ads_site)

        # Attach calculator
        slab.calc = calculator

        # Optimize
        opt = LBFGS(slab)
        opt.run(fmax=fmax, steps=steps)

        # Save final energy
        energy = slab.get_potential_energy()

        print(f"Final energy for {metal}: {energy:.6f} eV")
        f.write(f"{metal}\t{energy:.6f}\n")

print(f"\nDone. Energies saved to {outfile}")
