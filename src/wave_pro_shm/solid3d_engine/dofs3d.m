function iout = dofs3d(connectivity, iElement, nDofs, nElementNodes)
%
% Compute selection matrix used for assembly of iElement
% considering 3 degrees of freedom per node
%
% USAGE::
%
%   iout = dofs3d(connectivity,iElement,nDofs,nElementNodes)
%
% Arguments:
%     connectivity (integer):
%       connectivity matrix - list of degrees of freedom numbers in each element
%
%     iElement (integer):
%       element index
%
%     nDofs (integer):
%       number of degrees of freedom in the mesh
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
iout = logical(sparse(nElementNodes * 3, nDofs));
for i = 1:nElementNodes * 3
    iout(i, connectivity(iElement, i)) = 1;
end
% ---------------------------------------------------------------------------------------------------

end
