# Mixtures of Neural Network Experts with an Application to Phytoplankton Flow Cytometry Data

This is the code and data repository which supplements the article titled above. This method relies heavily upon the `flowmix` package developed by [Sangwon Hyun](https://github.com/sangwon-hyun/flowmix) and Jacob Bien.

To reproduce the results in the paper on a high-performance computing (HPC) environment, go into the `.slurm` scripts and specify any `#SBATCH` directives you need in order to make the scripts work with your specific resources, install any required R libraries (listed below) then follow these steps: 

```bash
Rscript 0_data_prep/00_corr.R

Rscript 1_simulation/10_sim_04_gen.R
chmod +x 1_simulation/11_sim_04.slurm
sbatch 1_simulation/11_sim_04.slurm
Rscript 1_simulation/12_results.R

Rscript 2_application/20_nested_cv/200_make_folds.R
chmod +x 2_application/20_nested_cv/201_pcr_oos_seq.sh
sbatch 2_application/20_nested_cv/201_pcr_oos_seq.sh
Rscript 2_application/20_nested_cv/202_pcr_oos_results.R

chmod +x 2_application/21_all_data/210_pcr_seq.sh
sbatch 2_application/21_all_data/210_pcr_seq.sh
Rscript 2_application/21_all_data/211_pcr_results.R
```

## Some Considerations:

- This workflow assumes you are using the root (top-level) directory as the working directory. If you don't have access to HPC, then the R scripts called by SLURM (`sim_04.R`, `pcr_oos.R`, and `pcr.R`) need to be called directly with the proper command-line arguments replacing the arguments passed by the SLURM scripts. If working on HPC, we recommend cloning this repository directly into the HPC environment, running all scripts there, and then using `scp` to copy the results and plot files to your PC for easier access.
- In the R scripts, file path separators are adjusted to your operating system automatically, but the SLURM files may need adjustment if you're not using a Unix-like operating system (for example, if you're on Windows, you may need to replace `/` with `\`).
- Don't run any `seedtab`-related scripts if you want to reproduce our results, since this will overwrite the seedtabs we used and therefore yield different initializations of the expectation-maximization algorithm.
- Finally, the `data_prep` script must be run first, but the `simulation`, `nested cv`, and `all_data` scripts may be run in any order or parallelized (they may take awhile to run).

## R Library Dependencies

- `flowmix`
- `flowtrend`
- `gridExtra`
- `ggplot2`
- `ggpubr`
- `RColorBrewer`
- `dplyr`
- `tidyr`
- `tibble`
- `magrittr`
- `lubridate`
- `reshape2`
- `parallel`
- `parallelly`
- `RhpcBLASctl`

## TODO

- test code
- clean `211_pcr_results.R` more