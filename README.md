# LEA3D
Computer-Aided Molecular Design program

[![badgepython](https://forthebadge.com/images/badges/made-with-python.svg)](https://www.python.org/downloads/release/python-370/)  [![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://chemoinfo.ipmc.cnrs.fr/)

**LEA3D** is a de novo design software which allows to optimize the structure of molecules. It is based on the following publications:

 [LEA3D: A Computer-Aided Ligand Design for Structure-Based Drug Design](https://pubs.acs.org/doi/10.1021/jm0492296)

 [e-LEA3D: a computational-aided drug design web server](https://academic.oup.com/nar/article/38/suppl_2/W615/1099650)

 [LEA (Ligand by Evolutionary Algorithm): A Genetic Algorithm for the Automated Generation of Small Organic Molecules](https://link.springer.com/article/10.1023/A:1008108423895)

LEA3D (Ligand by Evolutionary Algorithm) is designed to create new molecules by using a library of molecular fragments (structures in 3D) and by determining best combinations of molecular fragments that fit user-defined physicochemical properties (also called constraint function or fitness function). LEA3D is based on a genetic algorithm that evolves the molecular structures generation after generation until the emergence of fitted molecules. Each molecule of each generation is evaluated thanks to a fitness function (constraints) which can be either molecular properties, an affinity prediction by a docking program…

![example](/images/LEA3D.png)

Figure 1. General flowchart for LEA3D genetic algorithm.
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
	
Retrieve and unzip LEA3D repository in your desired folder. See below for running the program. The decompressed directory containing executables is called LEA3D-main.


## Run LEA3D
Create a folder called "Project", copy the content of the folder “examples” in Project directory and launch conda.
Create a folder "VISU" in the folder "Project".

 	conda activate lea3d

**1. Use lea3d to design molecules using some molecular properties of the aspirin molecular structure**

In this example, the file ligand-aspirin.in defines parameters of the genetic algorithm (a population of 10 molecules that will evolve over 30 generations) and indicates the fitness function file to read (ligand-aspirin.func). The fitness function includes 4 properties to evaluate (number of atoms, molecular weight, fsp3 value and the presence of 2 chemical functions (ester and acid)).

Execute:

	perl ../LEA3D-main/MAIN ligand-aspirin.in

The file summary.out indicates the score of candidate molecules of the last generation ranked in descending order. The file edesign.sdf contains all generated molecules over the run. The file VISU/list.sdf contains the best candidate molecule of each generation. The file popopop.txt contains the encoded molecules of the last generation. The file fitmoy.dat allows to plot the maximum, minimum and average scores in function of the generation number. In addition, the file operator.out records the crossover and mutation operations, indicates the difference in score value and which lego is involved (if any). It allows to analyze the efficiency of each operator. Of note, at the end of the run, the list of the privileged legos that improve candidate molecules is indicated.

Plot the maximum, minimum and average scores in function of the generation number

	python ../LEA3D-main/utils/plot-scores.py fitmoy.dat

 ![example](/images/plot.png)
 
 To visualize the best candidate of each generation:
 
 	pymol VISU/list.sdf

 ![example](/images/best-candidate.png)

The best candidate of the run is shown in the above Figure: it is the last molecule of the file list.sdf


**2. Use lea3d to generate 3 molecules using a pre-defined combination of legos**

In this example, the objective is to build molecules that are already encoded without evaluation of any properties. The file list_mol_sulfapyridine-aspirin_venetoclax contains the encoding for three molecules (one per line) and the file ligand-aspirin.in indicates which library of fragment to use (here, the default SDF file called all.sdf from the folder LEGO).

Execute:

	perl ../LEA3D-main/MAIN -v ligand-aspirin.in list_mol_sulfapyridine-aspirin_venetoclax
 
 To visualize the three generated molecules:
 
 	pymol mol_1.sdf mol_2.sdf mol_3.sdf

 ![example](/images/three-molecules.png)
 

**3. Use lea3d to evaluate molecules using the fitness function**

In this example, the objective is to use the program to evaluate the fitness function of a set of molecules. The SDF file of molecules is given as input and the file ligand-aspirin.in indicates the fitness function to use for the evaluation (ligand-aspirin.func).

Execute:

	perl ../LEA3D-main/MAIN -e ligand-aspirin.in three-molecules.sdf

The output is written on the screen:

![example](/images/result-lea3d-evaluation.png)

As indicated on the output, molecule number 2 in the sdf file has a score of 100%. This was expected as the aspirin itself is the second molecule of the file three-molecules.sdf. The file summary.txt indicates the score of screened molecules ranked in descending order of score.


**A more detailed description of the package can be found in the manual in docs folder (Parameter setting, How to create fragments, how to customize the core version of LEA3D...)**: [Manual-LEA3D-core.pdf](https://github.com/douguet/LEA3D/blob/main/docs/Manual-LEA3D-core.pdf)

## Licenses
LEA3D code is released under [the 3-Clause BSD License](https://opensource.org/licenses/BSD-3-Clause)

## Copyright
Copyright (c) 2018-2021, CNRS, Inserm, Dominique Douguet, All rights reserved.

## Reference

[Douguet D., Munier-Lehmann H., Labesse G. and Pochet S., LEA3D: A Computer-Aided Ligand Design for Structure-Based Drug Design, J. Med. Chem., 2005, 48, 2457-2468](https://pubs.acs.org/doi/10.1021/jm0492296). doi: 10.1021/jm0492296

[Douguet D., e-LEA3D: a computational-aided drug design web server, Nucleic Acids Res., 2010, 38, Suppl:W615-21](https://academic.oup.com/nar/article/38/suppl_2/W615/1099650). doi: 10.1093/nar/gkq322

[Douguet D., Thoreau E. and Grassy G., LEA (Ligand by Evolutionary Algorithm): A Genetic Algorithm for the Automated Generation of Small Organic Molecules, J. Comput.-Aided  Mol. Design, 2000, 14, 449-466.](https://link.springer.com/article/10.1023/A:1008108423895). doi: 10.1023/A:1008108423895
   
Bibtex format :

    @article{doi:10.1021/jm0492296,
    author = {Douguet, Dominique and Munier-Lehmann, Hélène and Labesse, Gilles and Pochet, Sylvie},
    title = {LEA3D: A Computer-Aided Ligand Design for Structure-Based Drug Design},
    journal = {Journal of Medicinal Chemistry},
    volume = {48},
    number = {7},
    pages = {2457-2468},
    year = {2005},
    doi = {10.1021/jm0492296},
    url = {https://pubs.acs.org/doi/10.1021/jm0492296},
    }
    
    @article{10.1093/nar/gkq322,
    author = {Douguet, Dominique},
    title = "{e-LEA3D: a computational-aided drug design web server}",
    journal = {Nucleic Acids Research},
    volume = {38},
    number = {suppl_2},
    pages = {W615-W621},
    year = {2010},
    month = {05},
    issn = {0305-1048},
    doi = {10.1093/nar/gkq322},
    url = {https://academic.oup.com/nar/article/38/suppl_2/W615/1099650},
    }
    
    @article{10.1023/A:1008108423895,
    author = {Douguet, Dominique and Thoreau, Etienne and Grassy, Gerard},
    title = "{A genetic algorithm for the automated generation of small organic molecules: Drug design using an evolutionary algorithm}",
    journal = {Journal of Computer-Aided Molecular Design},
    volume = {14},
    number = {},
    pages = {449-466},
    year = {2000},
    doi = {10.1023/A:1008108423895},
    url = {https://doi.org/10.1023/A:1008108423895},
    }


