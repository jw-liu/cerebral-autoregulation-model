function Rmatrix_ext = build_condition_network(net, cal, c)
%   Constructs the full extended cardio-cerebral resistance matrix (Fig.1c):
%     Heart ---> Aorta (R_aorta) --R_up[i]--> RICA / RVA / LVA / LICA
%                                --R_arm-->  Upper body (ground)
%                                --R_body--> Trunk / lower body (ground)
%   Inputs:
%     net - network geometry struct (from initialize_model)
%     cal - model parameters (R_aorta, R_up,R_arm, R_body)
%     c   - condition index (1, 2, or 3)
%   Output:
%     Rmatrix_ext - [Nodenum_ext x Nodenum_ext] extended resistance matrix

Rmatrix_ext = zeros(net.Nodenum_ext);
Rmatrix_ext(1:net.Nodenum, 1:net.Nodenum) = net.Rmatrix_patho{c};
R_up_c = cal.R_up(:, c);
for i = 1:net.nStart
    Rmatrix_ext(net.node_aorta, net.startNodes(i)) = R_up_c(i);
    Rmatrix_ext(net.startNodes(i), net.node_aorta) = R_up_c(i);
end

% Heart -> Aorta
Rmatrix_ext(net.node_heart, net.node_aorta) = cal.R_aorta;
Rmatrix_ext(net.node_aorta, net.node_heart) = cal.R_aorta;
% Aorta -> Arm
Rmatrix_ext(net.node_aorta, net.node_arm) = cal.R_arm;
Rmatrix_ext(net.node_arm,   net.node_aorta) = cal.R_arm;
% Aorta -> Body
Rmatrix_ext(net.node_aorta, net.node_body) = cal.R_body;
Rmatrix_ext(net.node_body,  net.node_aorta) = cal.R_body;

end
