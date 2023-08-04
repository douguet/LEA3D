# LEA3D
Computer-Aided Molecular Design program

[![badgepython](https://forthebadge.com/images/badges/made-with-python.svg)](https://www.python.org/downloads/release/python-370/)  [![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://chemoinfo.ipmc.cnrs.fr/)

**LEA3D** is a de novo design software which allows to optimize the structure of molecules. It is based on the following publications: [LEA3D: A Computer-Aided Ligand Design for Structure-Based Drug Design](https://pubs.acs.org/doi/10.1021/jm0492296); [e-LEA3D: a computational-aided drug design web server](https://pubs.acs.org/doi/10.1021/jm0492296) [LEA (Ligand by Evolutionary Algorithm): A Genetic Algorithm for the Automated Generation of Small Organic Molecules](https://link.springer.com/article/10.1023/A:1008108423895)

LEA3D (Ligand by Evolutionary Algorithm) is designed to create new molecules by using a library of molecular fragments (structures in 3D) and by determining best combinations of molecular fragments that fit user-defined physicochemical properties (also called constraint function or fitness function). LEA3D is based on a genetic algorithm that evolves the molecular structures generation after generation until the emergence of fitted molecules. Each molecule of each generation is evaluated thanks to a fitness function (constraints) which can be either molecular properties, an affinity prediction by a docking program…

![example](/images/LEA3D.png)

Figure 1. General flowchart for LEA genetic algorithm.
An initial population of candidate solutions is generated, usually, by random process but an option allows to start with a pool of molecules. The fitness of each candidate is evaluated via a fitness function (or score), which takes as input a candidate solution and returns a numeric score. Selection criteria are applied to choose candidates based on their fitness score for breeding. Breeding functions, crossover and mutations (suppress, add, replace or permutate a fragment), are applied to produce new solutions that replace the parent solutions. The cycle (or generation g) continues until convergence criteria is met (usually, solutions are no more improved). 

**Here, we provide the core of LEA3D that could be used or modified by adapting new fitness functions and/or libraries of fragments.**

**Documentation**: [Manual-LEA3D-core.pdf](https://github.com/douguet/LEA3D/blob/main/docs/Manual-LEA3D-core.pdf)

**Website**: A webserver dedicated to drug design is available at https://chemoinfo.ipmc.cnrs.fr/LEA3D/index.html

## Requirements

LEA3D uses **Perl and Python3.7**. Of note, JavaScript is only used to create html pages to display results. Here, we also use the RDKit package to calculate molecular properties. 


## Virtual environment for python with conda (for Windows for example)

Install conda or Miniconda from [https://conda.io/miniconda.html](https://conda.io/miniconda.html)  
Launch Anaconda Prompt, then complete the installation:

	conda update conda
	conda create -n lea3d
	conda activate lea3d
	conda install python=3.7 numpy
 	conda install perl
  	conda install -c conda-forge rdkit
   	conda install matplotlib

(Optional) Additional packages for visualization with PyMOL:

  	conda install -c schrodinger pymol
 
 If an error regarding PyQT occurs then:
 
 	conda install -c anaconda pyqt
	
Retrieve and unzip LEA3D repository in your desired folder. See below for running the program **lea3d**. The directory containing executables is called lea3d-main.


## Run LEA3D
Create a folder "Project", copy examples from the folder "examples" and launch conda

 	conda activate lea3d

**1. Use lea3d to design molecules using some molecular properties of the aspirin molecular structure**

Here, the file ligand-aspirin.func is read (its name is indicated in file ligand-aspirin.in). It contains 4 properties to evaluate (number of atoms, molecular weight, fsp3 value and 2 chemical functions (ester + acid)).

	perl ../lea3d-main/MAIN ligand-aspirin.in
 
 To visualize the best candidate of each generation:
 
 	pymol VISU/list.sdf
**1. Use lea3d to generate 3 molecules using the combination of legos**
	
	perl ../lea3d-main/MAIN -v ligand-aspirin.in list_mol_sulfapyridine-aspirin_venetoclax
 
 To visualize the three generated molecules:
 
 	pymol mol_1.sdf mol_2.sdf mol3.sdf

**2. Use lea3d to design molecules using some aspirin molecular properties**

Here, the file ligand-aspirin.func is read (its name is indicated in file ligand-aspirin.in). It contains 4 properties to evaluate (number of atoms, molecular weight, fsp3 value and 2 chemical functions (ester + acid)).

	perl ../lea3d-main/MAIN ligand-aspirin.in
 
 To visualize the best candidate of each generation:
 
 	pymol VISU/list.sdf

  
## Licenses
LEA3D code is released under [the 3-Clause BSD License](https://opensource.org/licenses/BSD-3-Clause)

## Copyright
Copyright (c) 2018-2021, CNRS, Inserm, Université Côte d'Azur, Dominique Douguet, All rights reserved.

## Reference
[Douguet D. and Payan F., SenSaaS: Shape-based Alignment by Registration of Colored Point-based Surfaces, *Molecular Informatics*, **2020**, 8, 2000081](https://onlinelibrary.wiley.com/doi/full/10.1002/minf.202000081). doi: 10.1002/minf.202000081
   
Bibtex format :

	@article{10.1002/minf.202000081,
	author 		= {Douguet, Dominique and Payan, Frédéric},
	title 		= {sensaas: Shape-based Alignment by Registration of Colored Point-based Surfaces},
	journal 	= {Molecular Informatics},
	volume 		= {39},
	number 		= {8},
	pages 		= {2000081},
	keywords 	= {Shape-based alignment, molecular surfaces, point clouds, registration, molecular similarity},
	doi 		= {https://doi.org/10.1002/minf.202000081},
	url 		= {https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202000081},
	eprint 		= {https://onlinelibrary.wiley.com/doi/pdf/10.1002/minf.202000081},
	year 		= {2020}
	}
