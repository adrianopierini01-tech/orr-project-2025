# Computational ORR project
This project uses the OpenCatalysis *oc20* machine learning potentials, which is part of  of the UMA models.
A pretrained model can be used within the ASE environment through the *Fairchem* code. It can be installed via pip: 
<pre>
pip install farichem-core
</pre>
Users are required to register themselves on Hugging Face: https://huggingface.co/facebook/UMA

Example usage of the *oc20* model as a calculator in ASE:
```python
from fairchem.core import pretrained_mlip, FAIRChemCalculator

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
calc = FAIRChemCalculator(predictor, task_name="oc20")
```

Credits to [Jim Furness](https://github.com/JFurness1) for the [EnergyLeveller](https://github.com/JFurness1/EnergyLeveller) script.

## Project references
1. J. Phys. Chem. B 2004, 108, 46, 17886–17892 (https://doi.org/10.1021/jp047349j)
2. Chem. Rev. 2018, 118, 5, 2302–2312 (https://doi.org/10.1021/acs.chemrev.7b00488)
3. ACS Catal. 2021, 11, 10, 6059–6072 (https://doi.org/10.1021/acscatal.0c04525)
