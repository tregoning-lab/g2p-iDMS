# In silico deep mutational scan

Interactive heatmaps of predicted antibody binding across every amino acid substitution at the residues of interest in antibody-antigen complexes.

[View the heatmaps](https://vgerardi23.github.io/g2p-iDMS/)

## Using the site

Pick an antibody and antigen from the dropdowns, then:

- switch between **fold change** and **K<sub>D</sub>**
- **find a mutation** - `N487G` for one cell, `487` for a whole position
- **effect at least** - hide anything changing binding by less than 2-, 5-, 10-, 25- or 100-fold
- sort the table by fold change, K<sub>D</sub>, position or substitution
- hover a cell for exact values, drag to zoom, save a PNG, download the CSV

Colour is a blue-white-red diverging scale on a log axis, centred on wild type: **blue is tighter binding, red is weaker**, and white is no change. Grey cells are the wild-type residue. Each scan's range is symmetric about wild type, so the same colour means the same change in every scan.

## Available data

| Antibody | Antigen | Positions | Mutations | Predicted wild-type K<sub>D</sub> |
|---|---|---|---|---|
| Nirsevimab | RSV-A F | 27 | 513 | 2.0 nM |
| Nirsevimab | RSV-B F | 27 | 513 | 1.5 nM |
| RB1 (Clesrovimab) | RSV F | 23 | 437 | 3.3 nM |
| Casirivimab | SARS-CoV-2 Delta RBD | 23 | 437 | 6.5 nM |

## References
Antigen protein structures were produced using [RosettaCM](https://docs.rosettacommons.org/docs/latest/application_documentation/structure_prediction/RosettaCM). Binding affinity (K<sub>D</sub>) was predicted with [PRODIGY](https://rascar.science.uu.nl/prodigy/).
****


