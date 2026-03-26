function Psolutions = Pressuresolution(Rmatrix, LineNodes, P)
%   For each time step, solves the linear system  B * P = S
%   where B is the conductance matrix for free (interior) nodes and
%   S contains contributions from fixed-pressure (boundary) nodes.
%   Uses LU decomposition (pre-computed once, reused for all time steps).
%   Inputs:
%     Rmatrix   - [N x N] symmetric resistance matrix between nodes
%     LineNodes - [L x 2] line segment endpoint indices
%     P         - struct with boundary conditions:
%                 P.start(i).num   = inlet node index
%                 P.start(i).value = [Nt x 1] pressure time series
%                 P.end(i).num     = outlet node index
%                 P.end(i).value   = [Nt x 1] pressure time series
%                 P.inter(i).num   = intermediate node index
%
%   Output:
%     Psolutions - [Nodenum x Nt] pressure at each node for each time step

Nodenum = size(Rmatrix, 1); Nt = length(P.start(1).value);
nStart = length(P.start); nEnd= length(P.end); nInter = length(P.inter);
fixedNodes  = zeros(nStart + nEnd, 1); fixedValues = zeros(nStart + nEnd, Nt);
for i = 1:nStart
    fixedNodes(i) = P.start(i).num;
    fixedValues(i, :) = P.start(i).value(:)';
end
for i = 1:nEnd
    fixedNodes(nStart + i)= P.end(i).num;
    fixedValues(nStart + i, :) = P.end(i).value(:)';
end

freeNodes = zeros(nInter, 1);
for i = 1:nInter
    freeNodes(i) = P.inter(i).num;
end
nFree  = length(freeNodes);  nFixed = length(fixedNodes);

% --- Build conductance matrix from Rmatrix ---
Gmatrix = zeros(Nodenum); Linenum = size(LineNodes, 1);
for i = 1:Linenum
    n1 = LineNodes(i, 1);  n2 = LineNodes(i, 2); R= Rmatrix(n1, n2);
    if R > 0
        G = 1 / R; Gmatrix(n1, n2) = G;  Gmatrix(n2, n1) = G;
    end
end
B = zeros(nFree);
for i = 1:nFree
    ni = freeNodes(i); B(i, i) = sum(Gmatrix(ni, :)); 
    for j = 1:nFree
        if j ~= i
            nj = freeNodes(j); B(i, j) = -Gmatrix(ni, nj);      
        end
    end
end

% Pre-compute LU decomposition (B is constant across time steps)
[L_mat, U_mat, Perm] = lu(B); G_fixed = zeros(nFree, nFixed);
for i= 1:nFree
    for k = 1:nFixed
        G_fixed(i, k) = Gmatrix(freeNodes(i), fixedNodes(k));
    end
end

% --- Solve for each time step ---
Psolutions = zeros(Nodenum, Nt);
for t = 1:Nt
    for k = 1:nFixed
        Psolutions(fixedNodes(k), t) = fixedValues(k, t);
    end
    S = G_fixed * fixedValues(:, t); P_free= U_mat \ (L_mat \ (Perm * S));
    for i = 1:nFree
        Psolutions(freeNodes(i), t) = P_free(i);
    end
end

end
