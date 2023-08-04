# LEA3D
Computer-Aided Molecular Design program

[![badgepython](https://forthebadge.com/images/badges/made-with-python.svg)](https://www.python.org/downloads/release/python-370/)  [![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://chemoinfo.ipmc.cnrs.fr/)

**LEA3D** is a de novo design software which allows to optimize the structure of molecules. It is based on the publications: [LEA3D: A Computer-Aided Ligand Design for Structure-Based Drug Design](https://pubs.acs.org/doi/10.1021/jm0492296); [e-LEA3D: a computational-aided drug design web server](https://pubs.acs.org/doi/10.1021/jm0492296) [LEA (Ligand by Evolutionary Algorithm): A Genetic Algorithm for the Automated Generation of Small Organic Molecules](https://link.springer.com/article/10.1023/A:1008108423895)

Here, we provide the core of LEA3D that could be used or modified.  

![example](/images/LEA3D.png)

**Documentation**: [Manual-LEA3D-core.pdf](https://github.com/LEA3D/lea3d/blob/main/docs/Manual-LEA3D-core.pdf)

**Website**: A webserver dedicated to drug design is available at https://chemoinfo.ipmc.cnrs.fr/LEA3D/index.html

## Requirements

LEA3D uses **Perl and Python3.7**. Of note, JavaScript is only used to create html pages to display results.


## Virtual environment for python with conda (for Windows for example)

Install conda or Miniconda from [https://conda.io/miniconda.html](https://conda.io/miniconda.html)  
Launch Anaconda Prompt, then complete the installation:

	conda update conda
	conda create -n sensaas
	conda activate sensaas
	conda install python=3.7 numpy
 
Once Open3D downloaded:
  
 	conda install open3d-0.12.0-py37_0.tar.bz2

(Optional) Additional packages for visualization with PyMOL:

  	conda install -c schrodinger pymol
 
 If an error regarding PyQT occurs then:
 
 	conda install -c anaconda pyqt
	
Retrieve and unzip SENSAAS-PY repository in your desired folder. See below for running the program **sensaas.py** or **meta-sensaas.py**. The directory containing executables is called sensaas-py-main.

## Linux

Install:

1. Python3.7 and numpy
2. Open3D version 0.12.0 (more information at [http://www.open3d.org/docs/release/getting_started.html](http://www.open3d.org/docs/release/getting_started.html))

(Optional) Install additional packages for visualization with PyMOL:

3. PyMOL (a molecular viewer; more information at [https://pymolwiki.org](https://pymolwiki.org))
  
Retrieve and unzip SENSAAS-PY repository. The directory containing executables is called sensaas-py-main.

## MacOS

	Not tested

## Run Sensaas
To align a Source molecule on a Target molecule, the syntax is:
	
	sensaas.py sdf molecule-target.sdf sdf molecule-source.sdf slog.txt optim
	
Example:

	sensaas.py sdf examples/IMATINIB.sdf sdf examples/IMATINIB_mv.sdf slog.txt optim
	
You may have to run the script as follows:

	python sensaas.py sdf examples/IMATINIB.sdf sdf examples/IMATINIB_mv.sdf slog.txt optim


Don't worry if you get the following warning from Open3D: "*Open3D WARNING KDTreeFlann::SetRawData Failed due to no data.*". It is observed with conda on 	windows.

## Visualization 

You can use any molecular viewer. For instance, you can use PyMOL if installed (see optional packages) to load the Target and the aligned Source(s):

after aligning IMATINIB_mv.sdf on IMATINIB.sdf using sensaas.py:

	pymol examples/IMATINIB.sdf Source_tran.sdf 

or after executing meta-sensaas.py with several molecules:

	pymol examples/IMATINIB.sdf bestsensaas.sdf catsensaas.sdf
	
or after the post-processing:

	pymol examples/IMATINIB.sdf ordered-catsensaas.sdf


or after executing meta-sensaas.py with the repeat option (State 1 is Target and State 2 is the aligned Source):
	
	pymol examples/VALSARTAN.sdf sensaas-1.sdf

## Licenses
SENSAAS code is released under [the 3-Clause BSD License](https://opensource.org/licenses/BSD-3-Clause)

## Copyright
Copyright (c) 2018-2021, CNRS, Inserm, Université Côte d'Azur, Dominique Douguet and Frédéric Payan, All rights reserved.

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
