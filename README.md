# pafiber


## Description-branch main

Uses a molecular theory approach to compute the electrostatic potential and ion distribution around a peptide amphiphile (pa) fiber as described in Nap et al Front. Chem. vol 10, 2022.
The aqueous solution is characterized by pH and monovalent NaCl and RbCl.
The theory predict the charged found on self-assemlbed peptide amphiphilies as function of pH and ion concentration. In particular we computed the charge of the peptide amphiphiles nanofibers with sequence
c_16-V_2A_2E_2. The theory considers chargeregulation of the carboxylic groups arising form the acid-base equilibrium of the glutamic acid residues as well as ion condensation  to the glutamic acid residues. Likewise Ca-brindging is feasible. The dielectric constant is position dependendent and includes a electostrostatics self-energy or Born-energy.
The aqueous solution is characterized by pH and monovalent, NaCl, RbCl, KCl or CsCl.
Mixture of NaCl with either RbCl, KCl or CsCl can also be considered.
Likewise solution containing CaCl2 can be added considered. 
The theory is based on a Molecular Theory developed to described 
end-tethered weak polyelectrolytes, see e.g.,
Nap et al J. Polym. Sci. Part B Polym. Phys vol 44, p 2638-2662, 2005, J. Chem. Phys. vol 140 p 024910, 2014 and Soft Matter vol 14 p2365, 2018
Here, the spatial distribution of pa fibes  and glutamics acid residue density is considered to be fixed and inputted into the theory. Its spatial distribution is taken for MD -simualtion as described in  Nap et al Front. Chem. 2022.

### Prerequisites

The program uses the kinsol solver form the library package sundials.
Sundials/kinsol, version 2.6

### Installing

Modify Makefile such that in contains the appropriate linker flags to the sundials library.
Linker flags can be obtained form the sundials installation.



### Running

The program has the following input files. 
* input file 'input.in'   
* Input file 'concen.in'  contains salt concentration 
* Input file  'xpa.in' contains volume fraction of pa.
* Input file 'rhoEps.in' contains Glu density of pa.

Examples input files can be found in folder input_example <br/>
The program uses an general configuration input file called 'input.in' that contain key words
as descibed in myio.f90. The pa volume fraction and Glu number denisty are stored in "xpa.in" and ""rhoEpa.in" The particular input correspond to C_16V_2A_2E_2 with a line density of 17.3 1/nm.
See SI of  Nap et al Frontiers  2022.
First line indicates the number of salt concentrations to be considered. Subsequent lines are the values of the salt concentrations. <br/>
The variable "runtype" in input.in control a loop of pH and type of salt. <<br/>

Example file input_Na.in uses NaCl salt. <br/>
Example  file input_K.in uses KCL salt.  <br/>
Note to run program with KCl or CsCl use logical switch_Rb_with_K .true. or switch_Rb_with_C .true. in in input file. input.in.
This switch the volume for Rb to K or Cs respectively. All output file still have Rb in name. <br/>
In example folder contains a python program that computes the ion excess based on input file generated with the Fortran program.  

## Built With

* [Sundials](https://computation.llnl.gov/projects/sundials/)


## Versioning

version 1.3 09-23-2026

## Authors

* **Rikkert J Nap**
