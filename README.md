# Mixtures of Neural Network Experts with an Application to Phytoplankton Flow Cytometry Data

This is the code and data repository which supplements the article titled above. This method relies heavily upon the `flowmix` package developed by [Sangwon Hyun](https://github.com/sangwon-hyun/flowmix) and Jacob Bien.

To reproduce the results in the paper, please follow the order specified by the numbers prefixed to the directory and file names. For example, 

1. `0_data_prep/00_corr.R`
2. `1_simulation/10_sim_04_gen.R`
3. `1_simulation/11_sim_04.slurm`
4. `1_simulation/12_results.R`
5. `2_application/20_nested_cv/200_make_folds.R`
6. ...
7. `2_application/21_all_data/211_pcr_results.R`

Only the numbered files require the user to run them directly if using a high-performance computing environment. If you don't have access to high-performance computing, then the R scripts called by SLURM (`sim_04.R`, `pcr_oos.R`, and `pcr.R`) need to be called directly with the proper command-line arguments replacing the arguments passed by the SLURM scripts.

All scripts expect the root (top-level) directory to be the working directory. In the R scripts, file path separators are adjusted to your operating system automatically, but the SLURM files may need adjustment if you're not using a Unix-like operating system (for example, if you're on Windows, you may need to replace `/` with `\`).

Don't run any `seedtab`-related scripts if you want to reproduce our results, since this will overwrite the seedtabs we used and therefore yield different initializations of the expectation-maximization algorithm.

## TODO

- test code