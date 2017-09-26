# pcis - Toolbox for computing Polyhedral Controlled-Invariant Sets

This is an ongoing effort to combine code for polyhedral synthesis methods (subsume the repositories `cps-inv`, `multi-cinv` and parts of `mkz`). Targeted list of capabilities:

 - Pre algorithms via robust counterpart (x) and intersection
 - Invariance algorithms with termination guarantees (x)
 - Measurable and non-measurable disturbance
 - Disturbance affecting A matrix
 - State-dependent disturbance bounds (in conflict with disturbed A matrix) (x)

## TODO list

 - Implement things marked with (x) above
 - Include tests

## Longer term objectives

 - Implementation of `minHRep` that scales with size of output as opposed to input
 - `Pre` with pre-defined template with complexity bound
 - Algorithms for merging non-convex polyhedral unions