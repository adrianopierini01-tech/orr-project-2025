from ase.optimize import LBFGS
from ase.io import read
from fairchem.core import pretrained_mlip, FAIRChemCalculator

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calculator = FAIRChemCalculator(predictor, task_name="oc20")

# Settings
all_molecules = ["h2o","h2","o2"]
outfile = "data.gas.txt"

# Save the list of final energies
with open(outfile, "w") as f:
    f.write("Molecule\tFinal_Energy_eV\n")

    # Iterate over all molecules
    for mol in all_molecules:
        print(f"\n=== Optimizing {mol} in gas-phase ===")
    
        # Create the molecule
        xyz_path = "./molecules/"+molecule+".xyz"
        ... # HINT # read molecule
    
        # Associate the calculator
        ... # HINT # define the molecule's "calc" method
    
        # Set up LBFGS dynamics object
        ... # HINT # call the LBFGS optimizer
        opt.run(0.01, 100)
    
        # Write the final energy
        energy = molecule.get_potential_energy()
        print(f"Final energy for {mol}: {energy:.6f} eV")
        f.write(f"{mol}\t{energy:.6f}\n")
 
