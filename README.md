# md-cookbook

###intro:

This is a small library of python scripts and modules that I use to start [Folding@home](http://folding.stanford.edu) simulations. All modeling is done with [PyRosetta](http://www.pyrosetta.org/) and [ParmEd](http://parmed.github.io/ParmEd/html/index.html), and all simulations are done with [OpenMM](http://openmm.org/).

Enjoy!

---

###recipes:

  + **seed_from_pdb**: 

  Use a PDB file to solvate, minimize, and equilibrate a simulation.
  
  + **seed_mutant_from_wt**: 
  
  Mutate a single amino acid in an existing PDB and start a round of equilibration.*

  + **seed_from_seq**:
  
  Solvate and equilibrate a protein simulation from a primary sequence.*



<p align="right">*requires PyRosetta and ParmEd</p>

---

###dependecies:

  + OpenMM
  + ParmEd (optional)
  + PyRosetta (optional)
