# demo_TMB_design

Rosetta XML scripts and python scripts developed for the de novo design of transmembrane beta-barrel nanopores as described in the followinf manuscript deposited to BiorXiv:

"Sculpting conducting nanopore size and shape through de novo protein design", 
Samuel Berhanu, Sagardip Majumder, Thomas Muntener, James Whitehouse, Carolin Berner, Asim K Bera, Alex Kang, Binyong Liang, G Nasir Khan, Banumathi Sankaran,  View ORCID ProfileLukas Tamm, David J Brockwel, Sebastian Hiller, Sheena Radford, David Baker,  View ORCID ProfileAnastassia Andreevna Vorobieva
doi: https://doi.org/10.1101/2023.12.20.572500

The provided demo input files were developed to design 12-strands nanopores with a shear number of 14 and variable shapes of beta-barrel cross-sections. The starting inputs (blueprints and constraints files used by Rosetta to assemble beta-barrel backbones are generated using the scripts provided in the generate_blueprint directory. The backbones are then assembled and sequences are designed, refined and screened using the following pipeline for each beta-barrel architecture (scripts and instructions are provided in the sub-directories):
1. assemble_backbones (assemble coarse-granined representations of the protein backbones)
2. best_500 (full-atom refinement of the backbones in preparation of the next steps)
3. round1_hbnet (design of the mortise-tenon folding motifs)
4. round2_core (first round of combinatorial design of the pore-exposed residues)
5. round2_surf (first round of combinatorial design of the lipid-exposed surface residues)
6. round3_core (second round of combinatorial design of the pore-exposed residues)
7. networks_analysis (filtering of the designs based on the polar networks stabilizing the mortise-tenon folding motifs)
8. round3_surf (second round of combinatorial design of the lipid-exposed surface residues)
9. round4_surf (last round of combinatorial design of the lipid-exposed surface residues to fine-tune the hydrophobicity)

The final models of the ordered designs are provided in the design_models subdirectory.

