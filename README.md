### In this repository the SteadierCom tool is put to the test in simulating two types of communities
Notebooks show results of analysis with SteadierCom of two case studies: 
1. gut microbiome,
2. _E. coli_ mutants.

### How to navigate this repository

Main information is in the 'notebooks' directory  (Notebooks 5.x.x).

#### Repository Structure

```text
SteadierCom-Case-Studies/
│
├── archive/          # Archived files and previous analyses
├── data/             # Experimental data and metabolic models
│   └── Guillen2021/       # Data related to Guillen et al. (2021) study (generated and original)
│   └── Venturelli2018/    # Data related to Venturelli et al. (2018) study (generated and original)
├── figures/          # Figures and visualizations
├── functions/        # Python functions
│   └── data/              # Functions for simulation results and data from the original studies
│   └── simulation/        # Functions for simulation and analysis using the models
│   └── analysis/          # Functions for analysing and illustrating simulation resutls
├── notebooks/        # Jupyter notebooks for analyses and simulations. Main figures are in notebooks 5.x.x
├── output/           # Generated results and output files
├── scripts_bash/     # Bash scripts for automated workflows, mostly for running SteadierCom
├── scripts_python/   # Python scripts for automated workflows, mostly for creating files for simulations
```

### An example of the results
**Illustration of SteadierCom simulation of mass exchange in a gut microbiome community**
![Analysis results](figures/circos/combined_NONE.png)
_Description_: (A) Bar plots and lineplot production/consumption and relative abundance for each member. From the top: First subplot shows relative abundance for each member, second shows total mass produced by each member, third shows the colors for each member, and fourth (bottom) shows total mass consumed by each member. Colors in the consumption/production bar plots represent the compound class. (B) Circos plot showing mass exchange between members. Colors around the circle match with the species. The arrows in between represent mass exchange, and are colored based on the compound class. 

**Publication status**: This repository contains code developed for ongoing research at NTNU. The associated manuscript is in preparation. Example figures are illustrative and do not necessarily represent final research results.
