function [elementCornerNodeIndexes, edge1NodeIndexes, edge2NodeIndexes, edge3NodeIndexes, ...
          edge4NodeIndexes] = element2d_numbering(nElementNodesX, nElementNodesY)
%
% Output corner node indexes and edge node indexes of 2D spectral element
%
% USAGE::
%
%   [elementCornerNodeIndexes, edge1NodeIndexes, edge2NodeIndexes, edge3NodeIndexes, ...
%    edge4NodeIndexes] = element2d_numbering(nElementNodesX, nElementNodesY)
%
% Arguments:
%     nElementNodesX (integer):
%       number of element nodes along local x axis (ksi)
%
%     nElementNodesY (integer):
%       number of element nodes along local y axis (eta)
%
% Returns:
%     elementCornerNodeIndexes (integer):
%       corner node indexes
%
%     edge1NodeIndexes (integer):
%       edge 1 node indexes
%
%     edge2NodeIndexes (integer):
%       edge 2 node indexes
%
%     edge3NodeIndexes (integer):
%       edge 3 node indexes
%
%     edge4NodeIndexes (integer):
%       edge 4 node indexes
%
%
% .. Note:: Zigzag numbering is used in 2D solid spectral elements
%    ::
%
%      local coordinates
%      -----------------
%      |      eta      |
%      |      ^        |
%      |      |        |
%      |      -->ksi   |
%      |               |
%      |               |
%      -----------------
%      corner nodes numbering
%      3---------------4
%      |               |
%      |               |
%      |               |
%      |               |
%      |               |
%      |               |
%      1---------------2
%      edge numbering
%      -------3->-------
%      |               |
%      |               |
%      ^               ^
%      4               2
%      |               |
%      |               |
%      -------1->-------
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------
elementCornerNodeIndexes = [1, nElementNodesX, nElementNodesX * nElementNodesY - nElementNodesX + 1, ...
                            nElementNodesX * nElementNodesY];
% edge nodes
edge1NodeIndexes = 2:nElementNodesX - 1;
edge2NodeIndexes = zeros(1, nElementNodesY - 2);
for k = 1:nElementNodesY - 2
    edge2NodeIndexes(k) = (k + 1) * nElementNodesX;
end
edge3NodeIndexes = nElementNodesX * nElementNodesY - nElementNodesX + 2:nElementNodesX * nElementNodesY - 1;
edge4NodeIndexes = zeros(1, nElementNodesY - 2);
for k = 1:nElementNodesY - 2
    edge4NodeIndexes(k) = k * nElementNodesX + 1;
end
% ---------------------------------------------------------------------------------------------------

end
