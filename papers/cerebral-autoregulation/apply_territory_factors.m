function Rmatrix_out = apply_territory_factors(Rmatrix, LineNodes, territory_lines, factors)
%   Inputs:
%     Rmatrix         - base resistance matrix
%     LineNodes       - [L x 2] endpoint indices
%     territory_lines - [6 x 1] line indices for the 6 terminal territories
%     factors         - [6 x 1] scaling factors
%   Output:
%     Rmatrix_out - modified resistance matrix

Rmatrix_out = Rmatrix;
for T = 1:length(territory_lines)
    ln = territory_lines(T);
    n1 = LineNodes(ln, 1); n2 = LineNodes(ln, 2);
    R_new = Rmatrix(n1, n2) * factors(T);
    Rmatrix_out(n1, n2)= R_new; Rmatrix_out(n2, n1) = R_new; 
end

end
