#!/bin/bash
# Run the WithCAM coupled solver in MATLAB (batch mode).
# Usage: ./run_all.sh
#
# Output:
#   results_WithCAM.mat       - converged results for 3 conditions
#   fig_WithCAM_baseline.png  - Baseline vs clinical data
#   fig_WithCAM_PCA.png       - PCA variant vs clinical data
#   fig_WithCAM_ACA.png       - ACA variant vs clinical data

set -e
cd "$(dirname "$0")"

if ! command -v matlab &> /dev/null; then
    echo "MATLAB not found in PATH"; exit 1
fi

echo "Running main_WithCAM_coupled_solver ..."
matlab -nodisplay -nosplash -nodesktop -r "try; main_WithCAM_coupled_solver; catch ME; disp(ME.message); disp(ME.stack); end; exit;"

echo "Done."
