function P = build_boundary_conditions(node_heart, endNodes_ext, interNodes_ext, Pinlet_Cstep, alpha)

Nt = size(Pinlet_Cstep, 1);
P = struct('start', [], 'end', [], 'inter', []);

% Pressure source: heart
P.start(1).num   = node_heart;
P.start(1).value = Pinlet_Cstep(:, 2) * alpha;

% Grounded outlets
for i = 1:length(endNodes_ext)
    P.end(i).num   = endNodes_ext(i);
    P.end(i).value = zeros(Nt, 1);
end

% Interior nodes (unknowns)
for i = 1:length(interNodes_ext)
    P.inter(i).num   = interNodes_ext(i);
    P.inter(i).value = zeros(Nt, 1);
end

end
