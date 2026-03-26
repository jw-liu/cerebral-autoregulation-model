function [Linedata, Nodedata, Nodenum, Linenum, NodeCL] = fun_LineNodeCL(Linedata, Nodedata)
%   Matches vessel endpoints to node positions using spatial coordinates.
%   Inputs/Outputs:
%     Linedata - vessel segment struct array (modified in place)
%     Nodedata - node struct array (modified in place)
%     Nodenum  - total number of nodes
%     Linenum  - total number of line segments
%     NodeCL   - [Nodenum x 4] node-to-line connectivity matrix

Nodenum= length(Nodedata); Linenum = length(Linedata);
NodeCL= zeros(Nodenum, 4);
Endspace1 = cell2mat(arrayfun(@(x) x.space(1,:),   Linedata, 'UniformOutput', false)');
Endspace2 = cell2mat(arrayfun(@(x) x.space(end,:), Linedata, 'UniformOutput', false)');

for i = 1:Nodenum
    cs = Nodedata(i).space;
    L1 = find(ismember(Endspace1, cs, 'rows'));
    L2 = find(ismember(Endspace2, cs, 'rows'));
    Nodedata(i).CL = [L1; L2]';
    NodeCL(i, 1:length(Nodedata(i).CL)) = Nodedata(i).CL;
    if ~isempty(L1)
        for j = 1:length(L1)
            Linedata(L1(j)).Endpoints(1) = i;
        end
    end
    if ~isempty(L2)
        for j = 1:length(L2)
            Linedata(L2(j)).Endpoints(2) = i;
        end
    end
end

end
