# TaskSwitching

# TaskSwitching — MATLAB analysis & plotting code

Analysis and plotting routines for the paper (release **v1.0.0**).  
Tested with **MATLAB R2024a**.

## Citation & DOIs
**Code (this release v1.0.0):** https://doi.org/10.5281/zenodo.17288627  
**Dataset:** https://doi.org/10.5281/zenodo.17288395

> For reproducibility, please cite the **version DOI** above.

## Overview
- Scripts to generate the main figures from the paper.
- Full dataset is archived on Zenodo (below).

## Requirements
- **MATLAB:** R2024a  
- **Toolboxes:** (none beyond base MATLAB for the provided plotting scripts; update here if you add dependencies)

## Quick start
```matlab
% 1) In MATLAB, open the repo root and add paths
addpath(genpath(pwd));

% 2) (Optional) Place the full dataset under data/raw (see "Data" below)

% 3) Reproduce figures
%   Exp1/Exp2 figures:
plot_Exp1Exp2(set "expID = 1");   % -> Fig.1F, Fig.2, Fig.3D,F,H–I, Supplementary 1–2
plot_Exp1Exp2(set "expID = 2");   % -> Fig.4F

%   Exp3 figures:
plot_Exp3();        % -> Fig.6

