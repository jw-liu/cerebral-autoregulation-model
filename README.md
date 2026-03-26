# WithCAM: Multiscale Cerebral Blood Flow Solver

MATLAB implementation of the CAM-incorporated 0D-1D coupled solver described in the manuscript "Multiscale modeling of blood circulation with cerebral autoregulation and network pathway analysis for hemodynamic redistribution in the vascular network with anatomical variations and stenosis conditions"

## Quick Start

```matlab
cd Multiscale_CAM_Code
main_WithCAM_coupled_solver   % one-click: runs all 3 conditions + generates figures
```

Or batch mode: `./run_all.sh`

## What It Does

Running `main_WithCAM_coupled_solver.m` performs:

1. **0D-1D coupled simulation with CAM** for three CoW configurations:
   - Condition 1: Baseline (complete CoW)
   - Condition 2: PCA (Fetal-type posterior cerebral artery, P1 segment absent)
   - Condition 3: ACA (Missing anterior cerebral artery A1 segment)

2. **WithoutCAM reference** computation via bifurcation propagation (without CAM)

3. **Comparison figure generation**: three bar charts comparing model results against clinical measurements (with/without CAM vs experimental mean ± std from Zarrinkoob et al. (2015) )

### Outputs

| File | Description |
| --- | --- |
| `results_WithCAM.mat` | Converged 0D flow, pressure, 1D area for all 3 conditions |
| `fig_WithCAM_baseline.png` | Baseline: WithCAM vs clinical data |
| `fig_WithCAM_PCA.png` | PCA variant: WithCAM & WithoutCAM vs clinical data |
| `fig_WithCAM_ACA.png` | ACA variant: WithCAM & WithoutCAM vs clinical data |

## Requirements

- MATLAB R2018b or newer versions (base MATLAB only; no additional toolboxes required) 
- No external data files beyond the three `.mat` input files included

## Code-to-Algorithm Mapping

| Algorithm step | Source file | Equation |
| --- | --- | --- |
| Initialization 
| Load geometry & parameters | `initialize_model.m` | Eq. (1) |
| Line-node connectivity | `fun_LineNodeCL.m` | -- |
| Extended network (Fig. 1c) | `build_condition_network.m` | Fig. 1c |
| Boundary conditions | `build_boundary_conditions.m` | -- |
| Outer loop: 0D-1D coupling 
| Outer loop driver | `main_WithCAM_coupled_solver.m` | -- |
| R update from 1D area | `update_resistance.m` | Eq. (9) |
| Outer convergence check | `compute_outer_convergence.m` | Eq. (10) |
| Middle loop: CAM autoregulation 
| CAM iteration | `solve_CAM_iteration.m` | Eqs. (5)-(6) |
| 0D pressure solver  | `Pressuresolution.m` | Eq. (2) |
| Flow from nodal pressure | `compute_Qdif.m` | Eq. (7) |
| Territory factor scaling | `apply_territory_factors.m` | -- |
| Comparison
| WithoutCAM | `compute_WithoutCAM.m` | -- |
| Relative flow computation | `compute_rel_flows.m` | -- |
| Comparison bar charts | `plot_condition_comparison.m` | -- |

## Data Files

### Input data

**`GeometryData.mat`** (20 KB)

| Variable | Size | Description |
| --- | --- | --- |
| `Pinlet` | 2001 x 2 | Inlet pressure waveform. Col 1: time ; Col 2: pressure |
| `startNodes` | 4 x 1 | Node indices of the four source inlets: RICA, RVA, LVA, LICA |
| `endNodes` | 45 x 1 | Node indices of the 45 terminal outlets across six territories |

**`ScanIPData.mat`** (732 KB)

Vessel geometry reconstructed from medical imaging provided by Ii et al. 2020.

| Field | Size | Description |
| --- | --- | --- |
| `Linedata` | struct (L entries) | Per-vessel: centreline (`CL`), `Endpoints`, `Radius`, `ScanIPname`, `space` |
| `Nodedata` | struct (N entries) | Per-node: connected line index (`CL`), `ScanIPname`, `space` |
| `Linenum` | scalar | Number of vessel segments (L = 130) |
| `Nodenum` | scalar | Number of nodes (N = 118) |
| `NodeCL` | 118 x 4 | Node-to-line connectivity matrix |
| `startNodes` | 1 x 4 | Source inlet node indices |
| `endNodes` | 1 x 45 | Terminal outlet node indices |

**`model_params.mat`** (1.1 MB)

Cardio-cerebral coupling parameters, CAM autoregulation parameters.

| Field | Size | Description |
| --- | --- | --- |
| `R_arm` | scalar | Aorta to upper body resistance (Fig. 1c) |
| `R_body` | scalar | Aorta to trunk/lower body resistance (Fig. 1c) |
| `R_up` | scalar  | Upstream resistances  |
| `upstream_nodes` | 6 x 1 | Node indices of the 6 territory upstream nodes |
| `q_bar` | 6 x 1 | Baseline territory flow set-points (Eq. 5) |
| `P1_bar` | 6 x 1 | Baseline distal pressure set-points (Eq. 5) |
| `Rsa_bar` | 6 x 1 | Baseline small artery resistance |
| `Rv` | 6 x 1 | Venous resistance per territory |
| `Gq` | 6 x 3 | Autoregulation gain per territory and condition |
| `factor` | 6 x 3 | Converged CAM autoregulation factors |
| `Q_1D` | {3 x 1} cell | 1D flow per vessel per condition |
| `A_1D` | {3 x 1} cell | 1D cross-sectional area per condition |
| `dx_cells` | 130 x 1 | Cell spacing (cm) per vessel |
| `NCELL` | scalar | Number of spatial cells per vessel (8) |

### Output data

**`results_WithCAM.mat`** (2.0 MB)

Converged results. Contains struct `results` and `cond_names`.

| Field | Size per condition | Description |
| --- | --- | --- |
| `Qdif` | 130 x 201 | 0D vessel flow over one cardiac cycle |
| `Psol` | 118 x 201 | Nodal pressure over one cardiac cycle |
| `Q_1D` | 130 x 201 | 1D vessel flow (spatial average) |
| `A_1D` | 130 x 201 | 1D cross-sectional area (spatial average) |
| `cam_rel` | 3 x 10 | WithCAM relative flow at 10 key vessels (%) |
| `nocam_rel` | 3 x 10 | WithoutCAM relative flow at 10 key vessels (%) |
| `Rmatrix_conv` | {3 x 1} | Converged resistance matrix per condition |
