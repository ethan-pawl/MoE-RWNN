# Mixtures of Neural Network Experts with an Application to Phytoplankton Flow Cytometry Data

This is the code and data repository which supplements the paper titled above.

## TODO

- gather seedtabs from HB
- confirm seedtabs are correct
- clean up code
  - simulation
    - GitHub doesn't track empty folders, so need to go back and call `dir.create`
    - `results.R`
    - replace all references `res` to the simulation data with `simdata`
    - replace all references `ints` to the intercepts with `clust2_intercepts`
    - get seedtabs from HB, confirm that they match
  - application
- rename half to nl
- add code where I create the nested CV folds