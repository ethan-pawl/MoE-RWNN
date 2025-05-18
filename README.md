# Mixtures of Neural Network Experts with an Application to Phytoplankton Flow Cytometry Data

This is the code and data repository which supplements the paper titled above.

## TODO

- gather seedtabs from HB
- confirm seedtabs are correct
- keep in mind GitHub doesn't track empty folders, so may need to call `dir.create` at some times
- clean up code
  - simulation
      - or just put README files in those folders
      - ask Justin if I can put the paper data on GitHub
    - `results.R`
    - replace all references `res` to the simulation data with `simdata`
    - replace all references `ints` to the intercepts with `clust2_intercepts`
    - get seedtabs from HB, confirm that they match
  - application
- rename half to nl
- add code where I create the nested CV folds