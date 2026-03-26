function Qdif = compute_Qdif(Psol, Linenum, LineNodes, Rmatrix)
%   Q_ij(t) = [ P_i(t) - P_j(t) ] / R_ij
%   Inputs:
%     Psol      - [Nodenum x Nt] node pressures at each time step
%     Linenum   - number of line segments
%     LineNodes - [Linenum x 2] endpoint node indices for each segment
%     Rmatrix   - [Nodenum x Nodenum] resistance matrix
%   Output:
%     Qdif - [Linenum x Nt] flow through each segment

Nt = size(Psol, 2); Qdif = zeros(Linenum, Nt);

for i = 1:Linenum
    n1 = LineNodes(i, 1); n2 = LineNodes(i, 2);
    R= Rmatrix(n1, n2);
    if R > 0
        Qdif(i, :) = (Psol(n1, :) - Psol(n2, :)) / R;
    end
end

end
