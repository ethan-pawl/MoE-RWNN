# Mixtures of Neural Network Experts with Application to Phytoplankton Flow Cytometry Data

This is the code and data repository which supplements the article titled above. Implementation of this method was facilitated by the `flowmix` package developed by [Sangwon Hyun](https://github.com/sangwon-hyun/flowmix), Mattias Rolf Cape, François Ribalet, and Jacob Bien.

To reproduce the results in the paper on a high-performance computing (HPC) environment, go into the `.slurm` scripts and specify any `#SBATCH` directives you need in order to make the scripts work with your specific resources, install any required R libraries (listed below), then run the scripts in the order specified by the file name prefixes. For example, run

1. `Rscript 0_data_prep/00_EDA.R`
2. `Rscript 0_data_prep/01_data_prep.R`
3. `Rscript 1_simulation/10_sim_gen.R`
4. `sbatch 1_simulation/11_sim.slurm`
5. And so on

## Some Considerations:

- This workflow assumes you are using the root (top-level) directory as the working directory. If you don't have access to HPC, then the R scripts called by SLURM (`sim.R`, `pcr_oos.R`, and `pcr.R`) need to be called directly with the proper command-line arguments replacing the arguments passed by the SLURM scripts. If working on HPC, we recommend cloning this repository directly into the HPC environment, running all scripts there, and then using `scp` to copy the results and plot files to your PC for easier access.
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
- `lubridate`
- `reshape2`
- `parallel`
- `parallelly`
- `RhpcBLASctl`
- `gtools`
- `matrixStats`
- `mvtnorm`
- `remotes`
  - to install `flowmix` and `flowtrend`
- `flowmatch`
- `mgcv`
- `gifski`
- `magick`
- `ggrepel`
- `gganimate`
- `patchwork`
- `mixdistreg`
- `deepregression`
- `ellipse`
- `randomForestSRC`
