# Discontinuous Galerkin code to solve Regge-Wheeler in Reissner-Nordstrom spacetime

In addition, we use hyperboloidal layers, a special coordinate system to access null infinity
on our computational domain. The right boundary is mapped to null infinity.

To do this we use a compactification factor when changing the spatial coordinate and a height
function when changing the time coordinate.

## Description of folders:

1. **dG_layers_codes**:
 
   i. **setup**:         creates the entire DG setup. Used from the book.  
   ii. **custom_codes**: routines written for specific tasks/functionality by Manas.
   iii. **conversion**:  routines that calculate r coordinate from r* coordinate.
   iv. **potential**:    routine that calculate the potential term in the wave eqn.
   v. **ERN/RN case**:   solves wave eqn in the extremal and sub-extremal RN case.
   
2. **Stefanos_RN**: C and cuda code by Gaurav that he used to produce data. Important files are ```kern.cu```, ```main.c```, ```body.cu ```.

**Remark-1**: Layers can be attached to either the left (folder named ```/only_horizon```), right boundary or to both boundaries. This has been done in the folder ```/both_horizon_scri```. 


## Get Started

1. Run ```WaveDriver_RN.m``` inside ```/RN_case```. 

## Status

1. Tests for both cases are pending. Goal is to complete these on a later date.

## References

1.[An observational signature for extremal black holes](https://arxiv.org/abs/2307.03963)

2. [Toward exponentially-convergent simulations of extreme-mass-ratio inspirals: A time-domain solver for the scalar Teukolsky equation with singular source terms](https://arxiv.org/abs/2307.01349)

3. [Nodal Discontinuous Galerkin methods, Algorithms, Analysis and applications](https://link.springer.com/book/10.1007/978-0-387-72067-8)