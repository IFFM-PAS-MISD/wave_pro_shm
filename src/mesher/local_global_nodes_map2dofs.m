function [I_G_dofs, I_L_dofs] = local_global_nodes_map2dofs(I_G, I_L, isGPUavailable)
%
% Convert nodal numbers into degrees of freedom numbers (considering 3 dofs per node)
%
% USAGE::
%
%   [I_G_dofs,I_L_dofs] = local_global_nodes_map2dofs(I_G,I_L,isGPUavailable)
%
% Arguments:
%     I_G (integer):
%       matrix of global node indexes - corresponding to local nodes,
%       dimensions [m, nColumns]
%
%     I_L (integer):
%       matrix of local node indexes - corresponding to global nodes,
%       dimensions [m, nColumns]
%
%     isGPUavailable (logical):
%       flag for gpu availability, if true - computations are on gpu
%       if false - computations are on cpu
%
% Returns:
%     I_G_dofs (integer):
%       matrix of global dof indexes - corresponding to local dofs,
%       dimensions [3*m, nColumns]
%
%     I_L_dofs (integer):
%       matrix of local dof indexes - corresponding to global dofs,
%       dimensions [3*m, nColumns]
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

nColumns = size(I_G, 2);

I_G_dofs = zeros(3 * size(I_G, 1), nColumns);
I_G_dofs(1:3:end, :) = 3 * I_G - 2;
I_G_dofs(2:3:end, :) = 3 * I_G - 1;
I_G_dofs(3:3:end, :) = 3 * I_G;

I_L_dofs = zeros(3 * size(I_L, 1), nColumns);
I_L_dofs(1:3:end, :) = 3 * I_L - 2;
I_L_dofs(2:3:end, :) = 3 * I_L - 1;
I_L_dofs(3:3:end, :) = 3 * I_L;

if isGPUavailable
    I_G_dofs = gpuArray(I_G_dofs);
    I_L_dofs = gpuArray(I_L_dofs);
end

% ---------------------------------------------------------------------------------------------------

end
