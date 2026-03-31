function Rmatrix_out = apply_territory_factors(Rmatrix, LineNodes, territory_lines, factors)

Rmatrix_out = Rmatrix;
for T = 1:length(territory_lines)
    ln = territory_lines(T);
    n1 = LineNodes(ln, 1); n2 = LineNodes(ln, 2);
    R_new = Rmatrix(n1, n2) * factors(T);
    Rmatrix_out(n1, n2)= R_new; Rmatrix_out(n2, n1) = R_new; 
end

end
