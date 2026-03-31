function [Rmatrix_new, A_1D] = update_resistance(Qdif_conv, Q_1D, Rmatrix_current, net, cal, params, c)

Nt = size(Qdif_conv, 2); A_1D = cal.A_1D{c}; A_cells = cal.A_cells{c}; dx_cells = cal.dx_cells;
assert(size(Q_1D, 2) == Nt, '1D/0D time step mismatch');
assert(size(A_1D, 2) == Nt, '1D area time step mismatch');
C_R = params.C_R;  
R_eff = zeros(net.Nodenum_ext);
for iv = 1:net.Linenum
    Nodes = sort(net.Linedata(iv).Endpoints);
    dx = dx_cells(iv);   A_k = A_cells{iv};   
    Nt = size(A_k, 2); R_t = zeros(1, Nt);
    for nt = 1:Nt
        R_t(nt) = sum(dx ./ A_k(:, nt).^2) * pi^2 * C_R;
    end
    R_val = mean(R_t);
    R_eff(Nodes(1), Nodes(2)) = R_val;  R_eff(Nodes(2), Nodes(1)) = R_val;
end

Rmatrix_new = Rmatrix_current;
for iv = 1:net.Linenum
    ns = sort(net.Linedata(iv).Endpoints);
    R_old = Rmatrix_current(ns(1), ns(2));
    R_upd = (1 - params.gamma_relax)* R_old + params.gamma_relax * R_eff(ns(1), ns(2));
    Rmatrix_new(ns(1), ns(2)) = R_upd; Rmatrix_new(ns(2), ns(1)) = R_upd;
end

end
