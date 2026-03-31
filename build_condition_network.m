function Rmatrix_ext = build_condition_network(net, cal, c)

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
