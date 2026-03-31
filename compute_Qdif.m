function Qdif = compute_Qdif(Psol, Linenum, LineNodes, Rmatrix)

Nt = size(Psol, 2); Qdif = zeros(Linenum, Nt);
for i = 1:Linenum
    n1 = LineNodes(i, 1); n2 = LineNodes(i, 2);
    R= Rmatrix(n1, n2);
    if R > 0
        Qdif(i, :) = (Psol(n1, :) - Psol(n2, :)) / R;
    end
end

end
