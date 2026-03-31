function [Qdif_conv, Psol_conv] = solve_CAM_iteration(Rmatrix_current, net, cal, params, c)

nT   = 6;   % number of outlets originated into six territories, which refers to McConnell and Payne (2017)    
Gq_c = cal.Gq(:, c); territory_names = {'RACA','LACA','RMCA','LMCA','RPCA','LPCA'};

Vsa = ones(nT, 1); factor_pyn = ones(nT, 1); 
P1_prev = cal.P1_bar; q_hist= zeros(nT, nT + 1); 
Rmatrix_k = apply_territory_factors(Rmatrix_current, net.LineNodes_ext,net.territory_lines, factor_pyn);
Psol_k = Pressuresolution(Rmatrix_k, net.LineNodes_ext, net.P_bc);
Qdif_k = compute_Qdif(Psol_k, net.Linenum_ext, net.LineNodes_ext, Rmatrix_k);

for k = 1:params.max_cam
    for m = 1:nT
        q_m = mean(abs(Qdif_k(net.territory_lines(m), :)));
        dq  = (q_m - cal.q_bar(m)) / (cal.q_bar(m) + eps);
        if q_m <= cal.q_bar(m)
            dCa = cal.dCa_plus(m);       
        else
            dCa = cal.dCa_minus(m);     
        end
        Ca_m = cal.Ca_bar(m) + 0.5 * dCa * tanh(-2 * Gq_c(m) / (dCa + eps) * dq);
        P1_m = mean(abs(Psol_k(cal.upstream_nodes(m), :)));  dP1  = P1_m - P1_prev(m);
        Vsa(m) = Vsa(m) + Ca_m * dP1 / (cal.Rsa_bar(m) + eps);
        Vsa(m) = max(min(Vsa(m), 5.0), 0.2); P1_prev(m) = P1_m;
        factor_new = 1 / (Vsa(m)^2);
        factor_pyn(m) = factor_pyn(m) + params.relax_payne * (factor_new - factor_pyn(m));
        factor_pyn(m) = max(min(factor_pyn(m), 20.0), 0.05);
        Rmatrix_k = apply_territory_factors(Rmatrix_current, net.LineNodes_ext,net.territory_lines, factor_pyn);
        Psol_k= Pressuresolution(Rmatrix_k, net.LineNodes_ext, net.P_bc);
        Qdif_k = compute_Qdif(Psol_k, net.Linenum_ext, net.LineNodes_ext, Rmatrix_k);
    end
    q_now = zeros(nT, 1);
    for T = 1:nT
        q_now(T) = mean(abs(Qdif_k(net.territory_lines(T), :)));
    end
    col = mod(k - 1, nT + 1) + 1;
    if k > nT + 1
        dq_vec = abs(q_now - q_hist(:, col)); max_dq = max(dq_vec);
        thr = params.eps_cam * mean(cal.q_bar);
        if max_dq < thr
            fprintf('    CAM converged at k = %d  (max|dq| = %.2e)\n', k, max_dq);
            break;
        elseif k <= 20 || mod(k, 10) == 0
            fprintf('      k=%3d | max|dq| = %.4e  (threshold: %.2e)\n', k, max_dq, thr);
        end
    end
    q_hist(:, col) = q_now;
end

Rmatrix_final = apply_territory_factors(Rmatrix_current, net.LineNodes_ext,net.territory_lines, factor_pyn);
Psol_conv= Pressuresolution(Rmatrix_final, net.LineNodes_ext, net.P_bc);
Qdif_conv= compute_Qdif(Psol_conv, net.Linenum_ext, net.LineNodes_ext, Rmatrix_final);

end
