# Justification behind chosen competitors

- Reviewer explicitly mentioned GAMs, LASSO, Random Forest, and Deep Learning
  - GAM incorporates dimension reduction
  - Rugamer is probabilistic, is extremely similar in statistical framework, and incorporates deep learning in a statistically principled way
- MoE choice
  - Rugamer
    - Mentioned above
    - Also makes it clear that the change in the distribution is the target, not just the conditional mean (these two things are practically identical, but conceptually different)
  - Etienam
    - Very similar in modeling the conditional mean structure, but GP instead of NN
    - One of the most modern entries in the Tresp -> Rasmusses & Ghahramani -> Meeds & Osindero -> Yuan & Neubauer -> Nguyen & Bonilla -> Gadd -> ... line of literature
  - Liu
    - Does not model the distribution of the response data; is therefore not probabilistic
  - Gadd
    - Uses a DP (we assume the number of clusters is reasonably well-known *a priori*). This poses a theoretical challenge of how to merge clusters across time points (effective number of clusters could be very different across time points).
    - Claims scalability, but in practice is probably not scalable to our data size. 
    - Could incorporate dependence by putting y2 and y3 as x's since the model is joint, but y's and x's are fundamentally treated differently so this is not a good comparison.
    - Metropolis-Hastings and HMC steps challenge reproducibility.
    - 