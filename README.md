# TIS-DNA-Ex-Ion
OPENMM script to simulate coarse-grained DNA in explicit ion
1. Look at the Figure "TIS-DNA-construct.pdf" for a general TIS-DNA / TIS-ION model representation
2. To simulate a N nucleotide long ssDNA/dsDNA sequence, we add one dummy nucleotide at the 5' end and one dummy nucleotide at the 3' end of each chain 
3. The dummy nucleotide at the 5' end doesn't contain any phosphate group while the dummy nucleotide at the 3' end contains a phosphate group which is not charged
4. Sequence details of the dummy nucleotides do not matter as the bases do not take part in stacking/hydrogen bonding
5. In order to simulate 3 nucleotide long dA ssDNA, in the "user_input.py", specify the sequence with the two added dummy nucleotides: AAAAA. In the same file, mention if you are interested in simulating a ssDNA or dsDNA.
6. In the same "user_input.py" provide information regarding box-len and concentration of each type of ions. The present code lets you choose between Na, Mg, Ca. The number of ions are calculated based on the box length info and ion concentration. Make sure to check the final pdb before running the simulations.
7. On making changes to the file "user_input.py", run the script "generate_all_atm_cg.py"
8. Following that, run the packmol script to generate the final input pdb structure, "packmol<packmol_input.inp"
9. once the pdb file has been generated, run the command "python simulation_run.py"
10. go to Releases/v0.0.0-TIS-ION to find the updated code
