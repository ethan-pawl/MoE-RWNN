# Justification behind chosen competitors

- Reviewer explicitly mentioned GAMs, LASSO, Random Forest, and Deep Learning
  - GAM incorporates dimension reduction
  - Rugamer is probabilistic, is extremely similar in statistical framework, and incorporates deep learning in a statistically principled way
- MoE choice
  - Rugamer
    - Mentioned above
  - Liu
    - Does not model the distribution of the response data; is therefore not probabilistic
  - Gadd
    - Uses a DP (we assume the number of clusters is reasonably well-known *a priori*). This poses a theoretical challenge of how to merge clusters across time points (effective number of clusters could be very different across time points).
    - Claims scalability, but in practice is probably not scalable to our data size. 
    - Could incorporate dependence by putting y2 and y3 as x's since the model is joint, but y's and x's are fundamentally treated differently so this is not a good comparison.
    - Metropolis-Hastings and HMC steps challenge reproducibility.