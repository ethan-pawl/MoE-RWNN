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

# Spatial/Temporal Clustering

- [X] RandomNet
  - untrained CNN-LSTM -> k-means -> prune clustering -> ensemble
  - no covariates
  - not probabilistic
- [X] Discrimination/Clustering for MV Time Series
  - Reviwer seems very focused on multivariate setting
  - Review classification to know the difference b/w that and discrimination/clustering  
  - "we can be fairly sure of the correct partitions in advance"
    - assumes initial rough clustering or supervised scenario -> known spectral densities/categories or able to get an estimate
    - classifies an entire time series into a known label
    - no covariates, therefore the primary inferential objective is missing
    - FOR CLUSTERING:
      - don't need the labels in advance, just estimate the spectral densities and cluster based on pairwise divergences
      - but we still don't have a latent distribution for each cluster
- [X] SLEX (just one)
  - exact same problems as above paper
- [X] Model-Based Clustering of Multiple Time Series
  - true unsupervised problem
  - allows for covariates
  - univariate
- Multivariate Longitudinal
- [X] Spatial clustering of time series
  - Assume data are observed over an entire spatial domain at each time point (each location has a time series/each time point has a spatial data set)
  - No covariates, just plain AR
  - Constant probabilities over time
  - Assigns an entire time series a cluster membership so that regions of the spatial domain are clustered
  - Two-stage estimation
  - MRF assumes a lattice structure
- [X] Dynamic mixture of experts models for online prediction
  - Very close, but
    - Assumes same entities are observed at every time point
    - Assumes random-walk time-varying coefficients (may be too flexible)
      - We want to assume a bit more structure on the relationships
- [maybe] Classic switching regression reference
- [maybe] Variational learning for switching state-space models
  - Segments a time series (assumes states change across time, which doesn't make sense in our case even if the particles were shared across time)

- See about HMMs
  - Punzo - clarify what they're actually doing
- See spatial clustering overview
- Review functional LDA if reviewer mentioned it
- Also mention Rugamer and put everything in a new introduction paragraph
- Ignore point processes because that is definitely not our setting
- Put RandomNet-Rugamer (maybe Etienam) into ChatGPT and ask for a draft paragraph or two

# New Literature Review


