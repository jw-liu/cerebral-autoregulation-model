function rel = compute_rel_flows(Qdif, line_idx, idx_source, idx_sink)
%   Inputs:
%     Qdif       - [N x Nt] flow matrix (or [N x 1] vector)
%     line_idx   - [10 x 1] line indices for the 10 key vessels
%     idx_source - indices within line_idx for inlets (1:4)
%     idx_sink   - indices within line_idx for outlets (5:10)
%   Output:
%     rel - [1 x 10] relative flow percentage

Q = zeros(length(line_idx), 1);
for i = 1:length(line_idx)
    if line_idx(i) > 0
        if size(Qdif, 2) > 1
            Q(i) = abs(mean(Qdif(line_idx(i), :)));
        else
            Q(i) = abs(Qdif(line_idx(i)));
        end
    end
end

rel = zeros(1, length(line_idx));
s1 = sum(Q(idx_source)); s2 = sum(Q(idx_sink));
if s1 > 0, rel(idx_source) = Q(idx_source)' / s1 * 100; end
if s2 > 0, rel(idx_sink)   = Q(idx_sink)'   / s2 * 100; end

end
