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

2. **Comparison figure generation**: three bar charts comparing model results against clinical statistics (with/without CAM vs [mean ± s.d.] from Zarrinkoob et al. (2015) )

### Outputs

| File | Description |
| --- | --- |
| `fig_WithCAM_baseline.png` | Baseline: WithCAM vs clinical data |
| `fig_WithCAM_PCA.png` | PCA: WithCAM & WithoutCAM vs clinical data |
| `fig_WithCAM_ACA.png` | ACA: WithCAM & WithoutCAM vs clinical data |

## Requirements

- MATLAB R2018b or newer versions (base MATLAB only; no additional toolboxes required) 
- No external data files beyond the three `.mat` input files included

- ![demo](ArrowSt1Paths.gif)
