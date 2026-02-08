\# Reproducible R scripts for compost adoption analysis (EFA/EGA/CFA)



This repository contains the R scripts used to reproduce the analyses reported in our manuscript

(Compost practice example): polychoric correlations, parallel analysis, EGA + bootEGA, and CFA model comparisons.



\## DOI (archived release)

Zenodo DOI: 10.5281/zenodo.18517974



\## Requirements

\- R >= 4.2 (recommended: latest stable R)

\- Packages: psych, lavaan, EGAnet, qgraph, igraph, ggcorrplot, ggplot2, readxl



\## Data

Place the following files in `data/`:

\- `compost\_items.csv`  (item responses; columns must include the item names used in `R/00\_config.R`)

\- `cfa\_model\_syntax\_compost.xlsx` (Excel file containing 6-model lavaan syntax; sheet name in `R/00\_config.R`)



If the raw data cannot be shared, provide a `data/DATA\_README.md` describing:

\- variable names (items),

\- coding scheme (ordinal levels),

\- and how access can be requested or what derived objects are shared instead.



\## How to run (one command)

From the project root:

```r

source("run\_all.R")



