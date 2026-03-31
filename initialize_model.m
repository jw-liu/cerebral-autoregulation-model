function [net, cal, params] = initialize_model()

load('GeometryData.mat'); load('model_params.mat', 'cal');

[Linedata, Nodedata, Nodenum, Linenum, ~] = fun_LineNodeCL(Linedata, Nodedata);
LineNodes = vertcat(Linedata.Endpoints); nStart = length(startNodes);
Rmatrix = zeros(Nodenum);
for i = 1:Linenum
    Nodes = sort(Linedata(i).Endpoints);
    distances = sqrt(sum(diff(Linedata(i).space).^2, 2));  Len = sum(distances);
    Rmatrix(Nodes(1), Nodes(2)) = Len / mean(Linedata(i).Radius)^4 * cal.C_R;
end
Rmatrix1 = triu(Rmatrix, 1); Rmatrix = Rmatrix1 + Rmatrix1';  
vessel_names = cell(Linenum, 1);
for j = 1:Linenum
    if ~isempty(Linedata(j).ScanIPname) && ~isempty(Linedata(j).ScanIPname{1})
        vessel_names{j} = strtrim(Linedata(j).ScanIPname{1});
    else
        vessel_names{j} = sprintf('Line_%d', j);
    end
end

targetNames = {'RACA2','LACA2','RMCA2','LMCA2','RPCA2','LPCA2'}; territory_lines = zeros(6, 1);
for k = 1:6
    for j = 1:Linenum
        if strcmpi(vessel_names{j}, targetNames{k})
            territory_lines(k) = j; break;
        end
    end
end

labels_10 = {'RICA','LICA','RVA','LVA','RACA','LACA','RMCA','LMCA','RPCA','LPCA'};
idx_source = 1:4; idx_sink = 5:10;
line_idx_10 = zeros(10, 1);
for i = 1:10
    if i <= 4, name = labels_10{i}; else, name = [labels_10{i}, '2']; end
    for j = 1:Linenum
        if strcmpi(vessel_names{j}, name), line_idx_10(i) = j; break; end
    end
end

line_idx_10_cond = {line_idx_10, line_idx_10, line_idx_10};
line_idx_aca1 = line_idx_10;
line_idx_aca1(5) = find(cellfun(@(v) strcmpi(v,'RACA1'), vessel_names), 1);
line_idx_aca1(6) = find(cellfun(@(v) strcmpi(v,'LACA1'), vessel_names), 1);
line_idx_10_cond{3} = line_idx_aca1;

find_vn = @(name) find(cellfun(@(v) strcmpi(v, name), vessel_names), 1);
pv = struct('RICA',find_vn('RICA'),   'LICA',find_vn('LICA'), ...
    'RVA',find_vn('RVA'),'LVA',find_vn('LVA'),'BA',find_vn('BA'), ...
    'RMCA1',find_vn('RMCA1'), 'LMCA1',find_vn('LMCA1'), ...
    'RMCA2',find_vn('RMCA2'), 'LMCA2',find_vn('LMCA2'), ...
    'RACA1',find_vn('RACA1'), 'LACA1',find_vn('LACA1'), ...
    'RACA2',find_vn('RACA2'), 'LACA2',find_vn('LACA2'), ...
    'RPCA1',find_vn('RPCA1'), 'LPCA1',find_vn('LPCA1'), ...
    'RPCA2',find_vn('RPCA2'), 'LPCA2',find_vn('LPCA2'), ...
    'RPCOA',find_vn('RPCOA'), 'LPCOA',find_vn('LPCOA'), ...
    'ACOA',find_vn('ACOA'));

Rmatrix_A = Rmatrix;                                        % Baseline
pca_ln = 24;
n1p = LineNodes(pca_ln, 1); n2p = LineNodes(pca_ln, 2);
Rmatrix_B = Rmatrix;                                        % PCA
Rmatrix_B(n1p, n2p) = Rmatrix(n1p, n2p) * 1e8;
Rmatrix_B(n2p, n1p) = Rmatrix(n2p, n1p) * 1e8;
aca_ln = 21;
n1a = LineNodes(aca_ln, 1); n2a = LineNodes(aca_ln, 2);
Rmatrix_C = Rmatrix;                                        % ACA
Rmatrix_C(n1a, n2a) = Rmatrix(n1a, n2a) * 1e8;
Rmatrix_C(n2a, n1a) = Rmatrix(n2a, n1a) * 1e8;

node_aorta  = Nodenum + 1; node_heart  = Nodenum + 2;
node_arm = Nodenum + 3; node_body = Nodenum + 4;
Nodenum_ext = Nodenum + 4; Linenum_ext = Linenum + nStart + 3;

LineNodes_ext = [LineNodes; zeros(nStart + 3, 2)];
for i = 1:nStart
    LineNodes_ext(Linenum + i, :) = [node_aorta, startNodes(i)];
end
ln_source = Linenum + nStart + 1; ln_arm= Linenum + nStart + 2;
ln_body= Linenum + nStart + 3;
LineNodes_ext(ln_source, :) = [node_heart, node_aorta];
LineNodes_ext(ln_arm, :)= [node_aorta, node_arm];
LineNodes_ext(ln_body, :)  = [node_aorta, node_body];

startNodes_ext = node_heart; endNodes_ext   = [endNodes(:); node_arm; node_body];
interNodes_ext = setdiff(1:Nodenum_ext, [startNodes_ext; endNodes_ext]);

Pinlet_Cstep = Pinlet(1:10:end, :);
P_bc = build_boundary_conditions(node_heart, endNodes_ext, interNodes_ext,Pinlet_Cstep, cal.alpha);

nT = 6; territory_type = [3 3 1 1 2 2];  % ACA ACA MCA MCA PCA PCA
table_I = [4.71, 1.18, 0.08, 1.13, 0.06, 32.60, 12.18;   
    3.46, 0.87, 0.06, 0.83, 0.05, 60.50, 22.61;   
    3.96, 0.99, 0.07, 0.95, 0.05, 46.31, 17.30];  
Ca_bar = zeros(nT, 1); dCa_plus  = zeros(nT, 1); dCa_minus = zeros(nT, 1);
for T = 1:nT
    it= territory_type(T); Ca_bar(T)  = table_I(it, 3) /table_I(it, 1);
    dCa_plus(T)  = table_I(it, 4) / table_I(it, 1);
    dCa_minus(T) = table_I(it, 5) / table_I(it, 1);
end

mu_A = [36 36 14 15 11 12 21 21 8 8]; sig_A = [4 4 5 5 4 4 3 2 1 1];
mu_B = [40 37 10 13 10 13 21 19 8 8]; sig_B = [3 4 5 5 2 4 3 3 1 3];
mu_C = [27 45 14 14  0 19 22 22 7 7]; sig_C = [3 2 5 1 0 6 1 2 1 1];
norm_g = @(mu_v, sig_v, idx) deal(mu_v(idx)/sum(mu_v(idx))*100, sig_v(idx)/sum(mu_v(idx))*100);
exp_mu = zeros(3, 10); exp_sig = zeros(3, 10);
[a1,b1] = norm_g(mu_A, sig_A, idx_source); [a2,b2] = norm_g(mu_A, sig_A, idx_sink);
exp_mu(1,:) = [a1, a2]; exp_sig(1,:) = [b1, b2];
[a1,b1] = norm_g(mu_B, sig_B, idx_source); [a2,b2] = norm_g(mu_B, sig_B, idx_sink);
exp_mu(2,:) = [a1, a2]; exp_sig(2,:) = [b1, b2];
[a1,b1] = norm_g(mu_C, sig_C, idx_source); [a2,b2] = norm_g(mu_C, sig_C, idx_sink);
exp_mu(3,:) = [a1, a2]; exp_sig(3,:) = [b1, b2];

net = struct();
net.Linedata  = Linedata; net.Nodedata = Nodedata;
net.Nodenum= Nodenum; net.Linenum= Linenum;
net.LineNodes = LineNodes; net.startNodes= startNodes;
net.endNodes = endNodes; net.nStart = nStart;
net.vessel_names   = vessel_names; net.territory_lines = territory_lines;
net.line_idx_10  = line_idx_10; net.labels_10= labels_10;
net.idx_source = idx_source; net.idx_sink= idx_sink;
net.line_idx_10_cond = line_idx_10_cond; net.pv = pv;
net.Rmatrix_patho  = {Rmatrix_A, Rmatrix_B, Rmatrix_C};
net.Nodenum_ext = Nodenum_ext; net.Linenum_ext= Linenum_ext;
net.LineNodes_ext  = LineNodes_ext; net.node_aorta     = node_aorta;
net.node_heart= node_heart; net.node_arm  = node_arm;
net.node_body = node_body; net.P_bc = P_bc;

cal.Ca_bar= Ca_bar;
cal.dCa_plus= dCa_plus;
cal.dCa_minus = dCa_minus;

params = struct();
params.gamma_relax = 0.5;      
params.eps_cam= 5e-3;     
params.eps_outer   = 0.01;      
params.max_cam     = 500;       
params.max_outer   = 20;        
params.relax_payne = 0.3;      
params.C_R  = cal.C_R;  
params.exp_mu      = exp_mu;   
params.exp_sig     = exp_sig;   
params.cond_names  = {'Baseline', 'PCA (Fetal-type P1)', 'ACA (Missing or hypoplastic A1)'};
params.cond_titles = {'Baseline: Complete CoW', 'Fetal-type P1 segment of PCA','Missing or hypoplastic A1 segment of ACA'};
fprintf('[INIT] Done. Network: %d vessels, %d nodes\n', Linenum, Nodenum);
end
