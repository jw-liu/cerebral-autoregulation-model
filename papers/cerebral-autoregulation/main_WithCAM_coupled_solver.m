function  main_WithCAM_coupled_solver
%  CAM-incorporated 0D-1D multiscale coupled solver (as shown in Figure 2 in the manuscript).
%  Performs the iterative 0D-1D coupled simulation with cerebral
%  autoregulation for three CoW configurations:
%    Condition 1: Baseline (complete CoW)
%    Condition 2: PCA -- Fetal-type posterior cerebral artery (P1 absent)
%    Condition 3: ACA -- Missing anterior cerebral artery A1 segment
%  Outputs:
%    results_WithCAM.mat          -- converged 0D flow, pressure, 1D data
%    fig_WithCAM_baseline.png     -- Baseline comparison plot
%    fig_WithCAM_PCA.png          -- PCA comparison plot
%    fig_WithCAM_ACA.png          -- ACA comparison plot
%  Required data: GeometryData.mat, ScanIPData.mat, model_params.mat

clear; close all; clc; set(0, 'DefaultFigureRenderer', 'painters');
fprintf('\n==========================================================\n');
fprintf('  WithCAM Coupled 0D-1D Multiscale Solver\n');
fprintf('  Cerebral Artery Network with CoW Variants\n');
fprintf('==========================================================\n\n');

% STEP 0: Initialization
% NOTE on initial values and convergence speed:
% The baseline set-points in model_params.mat (q_bar, P1_bar, etc.) are taken from a prior converged simulation, so the initial state is already
% close to the solution.  This is why the CAM loop converges in ~8-10 iterations and the outer loop converges rapidly.  
% If different initial values are used (e.g. Q=0, P=0), more iterations will be needed.
% The algorithm supports up to max_outer=20 iterations with gamma=0.5 relaxation to ensure stable convergence from arbitrary starting points.
[net, cal, params] = initialize_model();
% Pre-allocate result storage
results= struct();
results.Qdif = cell(3, 1);   results.Psol= cell(3, 1);       
results.Q_1D = cell(3, 1);    results.A_1D = cell(3, 1);      
results.cam_rel = zeros(3, 10); results.nocam_rel = zeros(3, 10);
results.Rmatrix_conv = cell(3, 1);
cond_names = {'Baseline', 'PCA (Fetal-type P1)', 'ACA (Missing A1)'};

% CAM-incorporated 0D-1D Multiscale Model
for c = 1:3 % for three conditions
    fprintf('\n############################################################\n');
    fprintf('  CONDITION %d: %s\n', c, cond_names{c});
    fprintf('############################################################\n');
    Rmatrix_current = build_condition_network(net, cal, c);
    Q_prev = [];  A_prev = [];
    % 0D-1D coupling (cycle n), convergence: Eq.(10)
    % Each cycle: (1) solve 0D with CAM  ->  (2) solve 1D  ->  (3) update R
    % The outer loop iterates until Q and A converge between cycles.
    for n = 1:params.max_outer
        fprintf('\n  --- Outer cycle n = %d ---\n', n);
        [Qdif_conv, Psol_conv] = solve_CAM_iteration(Rmatrix_current, net, cal, params, c);
        Q_1D = cal.Q_1D{c}; dx_cells = cal.dx_cells;  NCELL = cal.NCELL;
        [Rmatrix_new, A_1D] = update_resistance(Qdif_conv, Q_1D, Rmatrix_current, net, cal, params, c);
        if n > 1 && ~isempty(Q_prev)
            J_n = compute_outer_convergence(Q_1D, Q_prev, A_1D, A_prev);
            fprintf('    J_%d = %.4e  (threshold: %.2e)\n', n, J_n, params.eps_outer);
            if J_n < params.eps_outer
                fprintf('    ==> Outer loop CONVERGED at n = %d\n', n);
                break;
            end
        end
        Q_prev = Q_1D;A_prev = A_1D; Rmatrix_current = Rmatrix_new; % update data
    end
    results.Qdif{c}= Qdif_conv(1:net.Linenum, :);
    results.Psol{c}  = Psol_conv(1:net.Nodenum, :);
    results.Q_1D{c}= Q_1D;  results.A_1D{c}= A_1D;
    results.Rmatrix_conv{c} = Rmatrix_current;
    results.cam_rel(c, :) = compute_rel_flows(results.Qdif{c}, net.line_idx_10_cond{c}, net.idx_source, net.idx_sink);
end
results.nocam_rel = compute_WithoutCAM(results, net);

% Generate comparison figures with simulated results for three conditions
fprintf('\n[*] Generating comparison figures ...\n');
fnames = {'fig_WithCAM_baseline', 'fig_WithCAM_PCA', 'fig_WithCAM_ACA'};
for c = 1:3
    plot_condition_comparison(results, net, params, c, fnames{c});
end

% Save results
save('results_WithCAM.mat', 'results', 'cond_names');
fprintf('\nResults saved to results_WithCAM.mat\n');
fprintf('\n==========================================================\n');
fprintf('  All done. Figures and results saved.\n');
fprintf('==========================================================\n');

end
