function rel = compute_rel_flows(Qdif, line_idx, idx_source, idx_sink)

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

rel = zeros(1, length(line_idx)); s1 = sum(Q(idx_source)); s2 = sum(Q(idx_sink));
if s1 > 0, rel(idx_source) = Q(idx_source)' / s1 * 100; end
if s2 > 0, rel(idx_sink)   = Q(idx_sink)'   / s2 * 100; end

end
