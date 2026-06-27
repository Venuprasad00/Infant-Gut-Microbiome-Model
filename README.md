## Overview

This repository contains the source code accompanying the manuscript:

> **A Mechanistic Model of Infant Gut Microbiome Development and its Modulation by Nutritional Interventions Following C-Section Delivery**

The model describes the temporal dynamics of four key microbial functional groups during the first six months of life while explicitly accounting for:

* Mode of delivery (vaginal vs. C-section)
* Oxygen depletion during early colonization
* Human milk oligosaccharides (HMOs)
* Galacto-oligosaccharides (GOS)
* Fructo-oligosaccharides (FOS)
* Dietary fiber introduced during weaning
* Microbial competition
* Nutritional interventions

The framework was developed to investigate microbiome restoration following C-section delivery.

---

## Repository Contents

| File                                 | Description                                                                                                                                            |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `ODEModelFile.jl`                    | Julia implementation of the mechanistic ordinary differential equation (ODE) model describing microbial population dynamics and substrate utilization. |
| `MLE_Parameter_Estimates_MATLAB.csv` | Maximum likelihood parameter estimates and 95% confidence intervals obtained from MATLAB parameter estimation.                                         |

---

## Model Components

The model tracks four microbial functional groups:

* **Bifidobacterium**
* **Pathogens** 
* **Bacteroides**
* **Clostridia**

together with five substrate pools:

* Human milk oligosaccharides (HMOs)
* Oxygen
* Dietary fiber
* Galacto-oligosaccharides (GOS)
* Fructo-oligosaccharides (FOS)

---

## Requirements

The model is implemented in **Julia** and requires:

* DifferentialEquations.jl
* CairoMakie.jl
* CSV.jl
* DataFrames.jl
* XLSX.jl

Additional MATLAB scripts were used for parameter estimation and confidence interval calculation.

---

## Citation

If you use this code in your work, please cite the accompanying manuscript.

