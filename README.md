# pcis - Toolbox for computing Polyhedral Controlled-Invariant Sets

This is an ongoing effort to combine code for polyhedral synthesis methods (subsume the repositories `cps-inv`, `multi-cinv` and parts of `mkz`). List of capabilities:

 - Pre algorithms via robust counterpart and intersection/projection
 - Iterative invariance algorithm with termination guarantee
   - State-dependent input bounds
   - Measurable and non-measurable disturbance
   - Parametrized/disturbed A matrix
   - State-dependent disturbance
 - One-shot search for implicit controlled-invariant sets
   - State-dependent input bounds

## Requirements
 - Recent version of Matlab and MPT 3.0 (http://control.ee.ethz.ch/~mpt).
 - Optional: Mosek for more efficient LP solving

To set up MPT to use Mosek:

``` 
 addpath /path/to/mosek/toolbox/r2014a
 savepath
 mpt_init 
```

## Installation
Add the `lib` folder to the Matlab path. Execute `runtests` in the `tests` folder to make sure everything works correctly. Optional: set up MPT to use a commercial solver.

## TODO list

 - Add 3D ACC example
 - Work out theory for w,v disturbance that is both p-dependent and x-dependent. Is there a conflict? 

## Longer term objectives

 - Integration of `minHRep` that scales better: https://github.com/mageecoe/indicate_nonredundant_halfplanes
 - `Pre` with pre-defined template with complexity bound
 - Algorithms for merging non-convex polyhedral unions
