# Mixtures of Neural Network Experts with Application to Phytoplankton Flow Cytometry Data

This is the code and data repository which reproduces the figures and numerical results in the paper "Mixtures of Neural Network Experts with Application to Phytoplankton Flow Cytometry Data", published in [*Environmetrics*](https://doi.org/10.1002/env.70148). Implementation of this method was facilitated by the `flowmix` package, developed by [Sangwon Hyun](https://github.com/sangwon-hyun/flowmix), Mattias Rolf Cape, François Ribalet, and Jacob Bien.

To reproduce the results in the paper on a high-performance computing (HPC) environment, go into the `.slurm` scripts and specify any `#SBATCH` directives you need in order to make the scripts work with your specific resources, install any required R libraries (listed below), then run the scripts in the order specified by the file name prefixes. For example, run

1. `Rscript 0_data_prep/0_EDA.R`
2. `Rscript 0_data_prep/1_data_prep.R`
3. `Rscript 1_simulation/0_sim_gen.R`
4. `sbatch 1_simulation/1_model_fitting.slurm`
5. And so on

## Some Considerations:

- This workflow assumes you are using the root (top-level) directory as the working directory. If you don't have access to HPC, then the R scripts called by SLURM (`1_simulation/model_fitting.R` and `2_application/3_PCs/cv.R`) need to be called directly with the proper command-line arguments, replacing the arguments passed by the SLURM scripts. If working on HPC, we recommend cloning this repository directly into the HPC environment, running all scripts there, and then using `scp` to copy the results and plot files to your PC for easier access.
- In the R scripts, file path separators are adjusted to your operating system automatically, but the SLURM files may need adjustment if you're not using a Unix-like operating system (for example, if you're on Windows, you may need to replace `/` with `\`). 

## R Library Dependencies

- `flowmix`
- `flowtrend`
- `gridExtra`
- `ggplot2`
- `ggpubr`
- `ggtext`
- `RColorBrewer`
- `maps`
- `dplyr`
- `tidyr`
- `tibble`
- `reshape2`
- `parallelly`
- `RhpcBLASctl`
- `matrixStats`
- `mvtnorm`
- `flowMatch`
- `mgcv`
- `mixdistreg`
- `deepregression`
- `randomForestSRC`
