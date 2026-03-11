# Multiple Stable States in Ecological Systems

Code repository accompanying a review paper on multiple stable states in ecological systems.

---

## Overview

This repository provides MATLAB scripts for simulating and visualising a two-species competition model that exhibits multiple stable equilibria. The example demonstrates how bistability arises from competitive exclusion dynamics, and produces publication-quality phase portraits and time-course plots.

---

## Model

The system is governed by a pair of coupled ODEs:

$$\frac{dN_1}{dt} = N_1(1 - N_1 - 2N_2) + \gamma_1$$

$$\frac{dN_2}{dt} = N_2(1 - 2N_1 - N_2) + \gamma_2$$

where $N_1$ and $N_2$ are population densities and $\gamma_1, \gamma_2 \geq 0$ are small immigration (or external input) terms that unfold the pitchfork bifurcation and create isolated stable branches.

---

## Scripts

| File | Description |
|---|---|
| `LV_example.m` | Main script: fixed-point analysis, time courses, phase portrait |

### What `LV_example.m` does

1. **Fixed-point analysis** — Locates all equilibria of the system using `fsolve` with a dense multi-start grid. Classifies each fixed point as *stable* (all eigenvalues of the Jacobian have negative real part) or *unstable/saddle* (at least one eigenvalue has positive real part). Results are printed to the command window.

2. **Figure 1 — Time courses** — Integrates the ODE from a single initial condition and plots $N_1(t)$ and $N_2(t)$.

3. **Figure 2 — Phase portrait** — Plots a normalised vector field (grey arrows) together with trajectories from a 7×7 grid of initial conditions (blue). Fixed points are overlaid:
   - 🔴 **Red filled circle** — stable fixed point
   - ⚫ **Black filled circle** — unstable fixed point / saddle

---

## Parameters

All adjustable parameters are declared at the top of `LV_example.m`:

| Parameter | Default | Description |
|---|---|---|
| `gamma1` | `0.1` | Immigration rate for species 1 |
| `gamma2` | `0.1` | Immigration rate for species 2 |
| `tspan` | `[0, 50]` | Integration time window |
| `IC_tc` | `[0.2, 0.1]` | Initial condition for the time-course plot |
| `nSide` | `7` | Grid size for phase-portrait initial conditions (`nSide²` trajectories) |
| `FS` | `22` | Base font size (pt) for all figures |
| `LW` | `2.0` | Trajectory line width |
| `FN` | `'Helvetica'` | Figure font (`'Times New Roman'` for serif) |

---

## Requirements

- MATLAB R2017b or later
- Optimization Toolbox (for `fsolve`)

No additional toolboxes are required.

---

## Usage

```matlab
% In MATLAB, navigate to the repository root and run:
LV_example
```

Fixed-point locations and their stability are printed to the command window. Two figure windows open automatically.

---

## License

MIT License — see [`LICENSE`](LICENSE) for details.
Copyright © 2026 Denis Patterson
