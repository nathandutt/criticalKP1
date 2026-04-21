# KP1 N-Soliton Solver

Generates N-soliton solutions to the KP1 equation using the τ-function expressed as a determinant, with pole dynamics accelerated via the Calogero-Moser system.

## Method

N-soliton solutions of KP1 are rational functions of x, y, t, so their structure is fully encoded by their poles. For a fixed t, the poles in x evolve as a function of y according to a **complex Calogero-Moser system** — an N-point ODE system. This avoids evaluating the full τ-function determinant on a dense grid and is significantly faster.

For each time step t, the program:

1. **Initialises pole positions** — diagonalises the Hirota matrix at `y_initial` to find the N pole positions in x.
2. **Evolves poles over y** — integrates the Calogero-Moser ODE system using RK4 across all y values.
3. **Reconstructs the field** — computes u(x, y) at all grid points from the pole positions via the rational ansatz.
4. **Locates critical points** — approximates critical points of u using bilinear interpolation on the reconstructed field.
5. **Stores results** — appends critical point data to a vector, written out as a `.npy` file.

## Input

Soliton parameters and initial conditions are read from `solitons.dat`.

## Output

Critical point positions (x, y) for each t, saved as a NumPy `.npy` file.

## Todo

1. Look into a faster solver for Calogero moser, maybe using hamiltonian/lax pair structure (Analog to symplectic solver)?
2. Maybe implement fast N body solver? (Fast multipole expansion/ Barnes Hut)
