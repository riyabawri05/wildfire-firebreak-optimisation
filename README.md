# Wildfire Firebreak Optimisation

**LSE MA324: Mathematical Modelling & Simulation, 2026**

## Problem

A fire authority manages a 21x21 grid landscape containing a settlement and wetland at risk from wildfires igniting in the southern scrubland. Their current plan places 17 firebreak cells in a horizontal east-west corridor but this was designed assuming fires travel northward, while observed winds blow predominantly from the west (70%) or south-east (30%). The goal is to find a better firebreak placement under a fixed budget.

## Approach

**Step 1: Network flow analysis.** I modelled fire propagation as a max-flow problem on the grid, with edge capacities set by spread probabilities. Shadow price analysis on the capacity constraints identifies the minimum cut, revealing the bottleneck corridors most critical to contain.

**Step 2: MILP optimisation.** I formulated two Mixed-Integer Linear Programming models in AMPL/Gurobi. Model A uses binary fire propagation, treating every active edge as carrying fire with certainty (worst-case). Model B uses continuous reachability, where fire probability decays multiplicatively along paths.

**Step 3: Monte Carlo simulation.** To evaluate robustness under uncertainty, I ran 10,000 scenarios sampling wind direction from a two-component von Mises mixture, wind speed from a Weibull distribution and ignition location uniformly across all vegetated cells, using acceptance-rejection sampling. Tail risk was evaluated using Conditional Value-at-Risk (CVaR) at the 90th percentile.

## Key finding

Worst-case scenarios are driven by ignitions inside the target zone, not by wind conditions. The recommended 13-cell plan matches the expected damage of the authority's 17-cell corridor at 29% better resource efficiency.

## Files

- `fire-simulator.r` — stochastic fire spread simulator provided by the module
- `landscape.csv` — 21x21 grid encoding vegetation and terrain
- `targets.csv` — settlement and wetland cell locations with damage weights
- `[FINAL] Cadmus Submission...` — full write-up including model derivations, code, and discussion

## Tools

R, AMPL(Gurobi)
