# MRP.py 
### A Parameterizer of Post-Translationally Modified Residues

This program streamlines the process of post-translationally modified residue parameterization.
MRP.py is designed to work with AmberTools and the Gaussian software package to handle molecular 
modeling and QM calculations, respectively. AmberTools15 or higher should be used in conjuction with this program. 

This code was developed by Patrick Sahrmann in the Goodwin research group at Auburn University. 

# Usage

MRP.py requires the anaconda and AMBER modules. The only initial atomic data required for MRP.py is the protein PDB file. MRP.py can be executed in the following manner:

`python MRP.py -i input_file.in -s step_number` 

Procedural details and input file construction information can be found in the corresponding paper. 

# License

This code is released under GNU General Public License V3.
