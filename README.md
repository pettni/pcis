# pcis - Toolbox for computing Polyhedral Controlled-Invariant Sets

This is an ongoing effort to combine code for polyhedral synthesis methods. Targeted list of capabilities:

 - Pre algorithms via robust counterpart and intersection
 - Invariance algorithms with termination guarantees
 - Measurable and non-measurable disturbance
 - Disturbance affecting A matrix
 - State-dependent disturbance bounds (in conflict with disturbed A matrix)


## Longer term objectives

 - Implementation of `minHRep` that scales with size of output as opposed to input
 - `Pre` with pre-defined template with complexity bound. 