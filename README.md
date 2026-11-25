# In silico deep mutational scan (*i*DMS)

Interactive heatmaps of predicted antibody binding across every amino acid substitution at the residues of interest in antibody-antigen complexes.

[View the heatmaps](https://tregoning-lab.github.io/g2p-iDMS/)

## Using the web page

Pick an antibody and antigen from the dropdowns, then:

- you can switch between **fold change** and **K<sub>D</sub>**
- you can **search a specific mutation** - eg. `N487G` for one cell, `487` for a whole position
- you can filter based on **effect** - hide anything changing binding by less than X-fold
- you can sort the table by fold change, K<sub>D</sub>, position or substitution
- you can hover a cell for exact K<sub>D</sub> and fold change values, drag to zoom, save a PNG, download the CSV

Colour is a blue-white-red diverging scale on a log axis, centred on wild type: **blue is tighter binding, red is weaker**, and white is no change. Grey cells are the wild-type residue. 

## Data available

| Monoclonal antibody | Antigen | Subtype or Strain | Positions mutated | Number of mutations | Predicted wild-type K<sub>D</sub> (nM) |
|---|---|---|---|---|---|
| Nirsevimab | RSV F | A | 27 | 513 | 2.0 |
| Nirsevimab | RSV F | B | 27 | 513 | 1.5 |
| RB1 (Clesrovimab) | RSV F | n/a | 23 | 437 | 3.3 |
| Casirivimab | SARS-CoV-2 RBD | Delta | 23 | 437 | 6.5 |

## How to cite

If you use this data or code, please cite:

> Mosscrop, L., Gerardi, V., Spearing, E., Kose, S., Talts, T., Thomas, C.,
> Rees-Spear, C., Paschos, K., Brown, J., Williams, T., Skinner, M.,
> de Lusignan, S., Zambon, M., Bravi, B. and Tregoning, J. (2026).
> An integrated in silico and in vitro genotype to phenotype pipeline to
> predict and characterise new RSV F site zero escape mutants. Preprint.
> doi:[10.21203/rs.3.rs-7884592/v1](https://doi.org/10.21203/rs.3.rs-7884592/v1)

## Licence

[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)
You may share and adapt the data, figures and code for any purpose, provided you give appropriate credit and cite the work above.


****


