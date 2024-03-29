function [I_G, I_L] = local_global_nodes_map(nElementNodes, elementNodes, nColumns)
%
% map local-global nodes and split into nColumns
% equivalent to mesh coloring algorithms but in respect to nodes instead of elements
% needed for gpu parallel computation
%
% USAGE::
%
%   [I_G,I_L] = local_global_nodes_map(nElementNodes,elementNodes,nColumns)
%
% Arguments:
%     nElementNodes (integer):
%       number of nodes in 3D solid spectral element
%
%     elementNodes (integer):
%       element nodes topology (element nodes connectivity)
%       dimensions [nElements in mesh, nNodes in element]
%
%     nColumns (integer):
%       number of columns into which node's map will be splitted
%
% Returns:
%     I_G (integer):
%       matrix of global node indexes - corresponding to local nodes
%
%     I_L (integer):
%       matrix of local node indexes - corresponding to global nodes
%
%
% (C) Copyright 2024 Piotr Fiborek, pfiborek@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

nNodes = max(max(elementNodes));
nElements = size(elementNodes, 1);
Psi = zeros(nNodes, 1);
phi_E = zeros(nNodes, nColumns);
phi_L = zeros(nNodes, nColumns);
for e = 1:nElements
    for i = 1:nElementNodes
        k = elementNodes(e, i);
        Psi(k) = Psi(k) + 1;
        phi_E(k, Psi(k)) = e;
        phi_L(k, Psi(k)) = i;
    end
end

I_G = zeros(nElementNodes * nElements / nColumns, nColumns);
I_L = zeros(nElementNodes * nElements / nColumns, nColumns);
H = zeros(nColumns, 1);

for k = 1:nNodes
    [~, S] = sort(H);
    for i = 1:Psi(k)
        H(S(i)) = H(S(i)) + 1;
        e = phi_E(k, i);
        I_G(H(S(i)), S(i)) = k;
        I_L(H(S(i)), S(i)) = (e - 1) * nElementNodes + phi_L(k, i);
    end
end

% ---------------------------------------------------------------------------------------------------

end
