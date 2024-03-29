function iout = dofs3dfifi(nodes, iElement, nNodes, nElementNodes)
%
% Compute selection matrix used for assembly of iElement
% considering 1 degree of freedom per node
%
% USAGE::
%
%   iout = dofs3d(nodes,iElement,nNodes,nElementNodes)
%
% Arguments:
%     nodes (integer):
%       connectivity matrix - list of nodes in each element
%
%     iElement (integer):
%       element index
%
%     nNodes (integer):
%       number of nodes in the mesh
%
%     nElementNodes (integer):
%       number of nodes in element
%
% Returns:
%     iout (logical):
%       sparse matrix with ones corresponding to degrees of freedom of iElement
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------
% Global degrees of freedom - selection matrix
iout = logical(sparse(nElementNodes, nNodes));
for i = 1:nElementNodes
    iout(i, nodes(iElement, i)) = 1;
end
% ---------------------------------------------------------------------------------------------------

end
