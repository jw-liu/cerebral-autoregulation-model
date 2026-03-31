function nocam_rel = compute_WithoutCAM(results, net)

pv = net.pv; Q_base = mean(abs(results.Qdif{1}), 2);
bv = struct(); fn = fieldnames(pv);
for k = 1:length(fn)
    bv.(fn{k}) = Q_base(pv.(fn{k}));
end
r_LPCA2 = bv.LPCA2 / (bv.LPCA2 + bv.LPCOA);
r_LPCOA = bv.LPCOA / (bv.LPCA2 + bv.LPCOA);
r_RMCA2 = bv.RMCA2 / (bv.RMCA2 + bv.RACA1);
r_RACA1 = bv.RACA1 / (bv.RMCA2 + bv.RACA1);
r_RACA2 = bv.RACA2 / (bv.RACA2 + bv.ACOA);
r_ACOA  = bv.ACOA  / (bv.RACA2 + bv.ACOA);
r_LMCA2 = bv.LMCA2 / (bv.LMCA2 + bv.LACA1);
r_LACA1 = bv.LACA1 / (bv.LMCA2 + bv.LACA1);

nocam_rel = zeros(3, 10);
nocam_rel(1, :) = results.cam_rel(1, :); q = Q_base;
q(pv.RPCA1) = 0;  q(pv.LPCA1) = bv.BA; q(pv.RPCA2) = 0; q(pv.RPCOA) = 0;
q(pv.LPCA2) = bv.BA * r_LPCA2;  q(pv.LPCOA) = bv.BA * r_LPCOA;
q(pv.RMCA1)= bv.RICA; q(pv.RMCA2) = bv.RICA * r_RMCA2;  q(pv.RACA1) = bv.RICA * r_RACA1;
q(pv.RACA2) =q(pv.RACA1) * r_RACA2;  q(pv.ACOA) = q(pv.RACA1) * r_ACOA;
q(pv.LMCA1) = bv.LICA + q(pv.LPCOA);
q(pv.LMCA2) = q(pv.LMCA1) * r_LMCA2;  q(pv.LACA1) = q(pv.LMCA1)* r_LACA1;
q(pv.LACA2) = q(pv.LACA1) + q(pv.ACOA);
nocam_rel(2,:) = compute_rel_flows(q, net.line_idx_10_cond{2}, net.idx_source, net.idx_sink);

q = Q_base;
q(pv.RACA1) = 0; q(pv.RMCA2) = bv.RMCA1; q(pv.RACA2) = 0;  q(pv.ACOA)  = 0;
q(pv.LACA2) = bv.LACA1; nocam_rel(3, :) = compute_rel_flows(q, net.line_idx_10_cond{3}, net.idx_source, net.idx_sink);

end
