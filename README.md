# Rapid Flash Calculation & Validation

This MATLAB project provides a high-performance implementation of flash calculations (Vapor-Liquid Equilibrium) for hydrocarbon mixtures, optimized for speed and thermodynamic consistency.

## Comparison vs. PVTi
The scripts benchmark the developed flash engine against industry-standard PVTi results, verifying:
- Phase fractions (Liquid/Vapor splits)
- Phase compositions ($x_i, y_i$)
- Equilibrium ratios ($K$-values)


## Requirements
- MATLAB (Optimization Toolbox recommended)
- PVTi output datasets (as .csv or .txt) for benchmarking
