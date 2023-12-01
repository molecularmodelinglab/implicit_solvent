# Re-ranking docked poses with AMBER + implicit solvent

We have used [GNINA](https://github.com/gnina/gnina) to propose poses for protein-ligand complex with a known crystal pose. It got a pose that's close to correct, but it didn't rank it first. Can our force field correctly rerank the poses?

We have a function that uses OpenMM to minimize a protein-ligand structure and return the minimized energy. We can also determine the Root-mean-squared-deviation (RMSD) between each docked pose and the known crystal pose. Low energy means the complex is more stable, and low RMSD means the complex is closer to the correct pose. Your task is to modify this code so that we can figure out:
1. Which pose that GNINA gave is the "correct" one? (lowest RMSD)
2. Is this pose also the pose with the lowest energy?
3. What does the closest docked pose look like compared to the crystal pose? Install [PyMol](https://pymol.org/2/) and use it to visualize the complex.

In order to get this code running, you'll need to create a new conda environment. Then you can install the required packages with:
```
conda install -c conda-forge -c omnia openmm openff-toolkit openff-forcefields openmmforcefields -y
pip install -r requirements.txt
```

All the starter code is in `implicit_solvent.py`. Good luck!
