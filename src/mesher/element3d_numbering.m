function [elementCornerNodeIndexes, Edge, Face, elementInteriorNodeIndexes] = ...
          element3d_numbering(nElementNodesX, nElementNodesY, nElementNodesZ)
%
% output node numbers in 8 corners, 12 edges, 6 faces and interior
% for complete high order solid spectral element
% number of element nodes along each axis can be different
%
% USAGE::
%
%   [elementCornerNodeIndexes,Edge,Face,elementInteriorNodeIndexes] = ...
%    element3d_numbering(nElementNodesX,nElementNodesY,nElementNodesZ)
%
% Arguments:
%     nElementNodesX (integer):
%       number of element nodes along local x axis (ksi)
%
%     nElementNodesY (integer):
%       number of element nodes along local y axis (eta)
%
%     nElementNodesZ (integer):
%       number of element nodes along local z axis (dzeta)
%
% Returns:
%     elementCornerNodeIndexes (integer):
%       corner node indexes
%
%     Edge (obj):
%       node indexes of each edge of spectral solid element
%
%     Face (obj):
%       indexes of each face of spectral solid element
%
%     elementInteriorNodeIndexes (integer):
%       interior node indexes
%
%
% .. Note:: Zigzag numbering is used in 3D solid spectral elements
%    ::
%
%      local coordinates
%                -----------------
%               / |             /|
%              /  |            / |
%             /   |           /  |
%            /    |          /   |
%           /     |         /    |
%          -----------------     |
%          |      |---^----|---- |
%          |    dzeta | eta|    /
%          |    /     |/   |   /
%          |   /      -->ksi  /
%          |  /            | /
%          | /             |/
%          -----------------
%      corner nodes numbering
%                7---------------8
%               / |             /|
%              /  |            / |
%             /   |           /  |
%            /    |          /   |
%           /     |         /    |
%          5---------------6     |
%          |      3--------|---- 4
%          |     /         |    /
%          |    /          |   /
%          |   /           |  /
%          |  /            | /
%          | /             |/
%          1---------------2
%     edge numbering
%                --------7--------
%               / |             /|
%              /  |            / |
%             8   |           6  |
%            /    12         /   11
%           /     |         /    |
%          --------5--------     |
%          |      |-----3--|---- |
%          |     /         |    /
%          |    /          |   /
%          9   4           10 2
%          |  /            | /
%          | /             |/
%          --------1--------
%     face numbering (vertical faces first, than bottom (5) and top (6))
%                -----------------
%               / |             /|
%              /  |            / |
%             /   |           /  |
%            /    |     4    /   |
%           /     |         /    |
%          -----------------  2  |
%          |   1  |--------|---- |
%          |     /         |    /
%          |    /  3       |   /
%          |   /           |  /
%          |  /            | /
%          | /             |/
%          -----------------
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

%% corner nodes
elementCornerNodeIndexes = [
                          1, nElementNodesX, nElementNodesX * nElementNodesY - nElementNodesX + 1, ...
                          nElementNodesX * nElementNodesY, ...
                          (nElementNodesZ - 1) * nElementNodesX * nElementNodesY + 1, ...
                          (nElementNodesZ - 1) * nElementNodesX * nElementNodesY + nElementNodesX, ...
                          (nElementNodesZ - 1) * nElementNodesX * nElementNodesY + nElementNodesX * nElementNodesY - nElementNodesX + 1, ...
                          nElementNodesZ * nElementNodesX * nElementNodesY];

%% edge nodes
% number of nodes at each edge
nEdgeNodes(1) = nElementNodesX - 2;
nEdgeNodes(2) = nElementNodesY - 2;
nEdgeNodes(9) = nElementNodesZ - 2;
nEdgeNodes(3) = nEdgeNodes(1);
nEdgeNodes(4) = nEdgeNodes(2);
nEdgeNodes(5) = nEdgeNodes(1);
nEdgeNodes(10) = nEdgeNodes(9);
nEdgeNodes(11) = nEdgeNodes(9);
nEdgeNodes(12) = nEdgeNodes(9);
nEdgeNodes(6) = nEdgeNodes(2);
nEdgeNodes(7) = nEdgeNodes(1);
nEdgeNodes(8) = nEdgeNodes(2);
% list of node numbers at each edge
Edge(8) = struct();
if nEdgeNodes(1) > 0
    Edge(1).elementEdgeNodeIndexes = 2:nElementNodesX - 1;
else
    Edge(1).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(2) > 0
    elementEdgeNodeIndexes = zeros(1, nEdgeNodes(2));
    for k = 1:nElementNodesY - 2
        elementEdgeNodeIndexes(k) = (k + 1) * nElementNodesX;
    end
    Edge(2).elementEdgeNodeIndexes = elementEdgeNodeIndexes;
else
    Edge(2).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(3) > 0
    Edge(3).elementEdgeNodeIndexes = nElementNodesX * nElementNodesY - nElementNodesX + 2:nElementNodesX * nElementNodesY - 1;
else
    Edge(3).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(4) > 0
    elementEdgeNodeIndexes = zeros(1, nEdgeNodes(4));
    for k = 1:nElementNodesY - 2
        elementEdgeNodeIndexes(k) = k * nElementNodesX + 1;
    end
    Edge(4).elementEdgeNodeIndexes = elementEdgeNodeIndexes;
else
    Edge(4).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(5) > 0
    Edge(5).elementEdgeNodeIndexes = nElementNodesX * nElementNodesY * (nElementNodesZ - 1) + 2: ...
                                     nElementNodesX * nElementNodesY * (nElementNodesZ - 1) + nElementNodesX - 1;
else
    Edge(5).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(6) > 0
    elementEdgeNodeIndexes = zeros(1, nEdgeNodes(6));
    for k = 2:nElementNodesY - 1
        elementEdgeNodeIndexes(k - 1) = nElementNodesX * nElementNodesY * (nElementNodesZ - 1) + nElementNodesX * k;
    end
    Edge(6).elementEdgeNodeIndexes = elementEdgeNodeIndexes;
else
    Edge(6).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(7) > 0
    Edge(7).elementEdgeNodeIndexes = nElementNodesX * nElementNodesY * nElementNodesZ - nElementNodesX + 2: ...
                                     nElementNodesX * nElementNodesY * nElementNodesZ - 1;
else
    Edge(7).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(8) > 0
    elementEdgeNodeIndexes = zeros(1, nEdgeNodes(8));
    for k = 1:nElementNodesY - 2
        elementEdgeNodeIndexes(k) = nElementNodesX * nElementNodesY * (nElementNodesZ - 1) + nElementNodesX * k + 1;
    end
    Edge(8).elementEdgeNodeIndexes = elementEdgeNodeIndexes;
else
    Edge(8).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(9) > 0
    elementEdgeNodeIndexes = zeros(1, nEdgeNodes(9));
    for k = 2:nElementNodesZ - 1
        elementEdgeNodeIndexes(k - 1) = (k - 1) * nElementNodesX * nElementNodesY + 1;
    end
    Edge(9).elementEdgeNodeIndexes = elementEdgeNodeIndexes;
else
    Edge(9).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(10) > 0
    elementEdgeNodeIndexes = zeros(1, nEdgeNodes(10));
    for k = 2:nElementNodesZ - 1
        elementEdgeNodeIndexes(k - 1) = (k - 1) * nElementNodesX * nElementNodesY + nElementNodesX;
    end
    Edge(10).elementEdgeNodeIndexes = elementEdgeNodeIndexes;
else
    Edge(10).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(11) > 0
    elementEdgeNodeIndexes = zeros(1, nEdgeNodes(11));
    for k = 2:nElementNodesZ - 1
        elementEdgeNodeIndexes(k - 1) = k * nElementNodesX * nElementNodesY;
    end
    Edge(11).elementEdgeNodeIndexes = elementEdgeNodeIndexes;
else
    Edge(11).elementEdgeNodeIndexes = [];
end
if nEdgeNodes(12) > 0
    elementEdgeNodeIndexes = zeros(1, nEdgeNodes(12));
    for k = 2:nElementNodesZ - 1
        elementEdgeNodeIndexes(k - 1) = k * nElementNodesX * nElementNodesY - nElementNodesX + 1;
    end
    Edge(12).elementEdgeNodeIndexes = elementEdgeNodeIndexes;
else
    Edge(11).elementEdgeNodeIndexes = [];
end
%% face nodes
Face(6) = struct();
if nEdgeNodes(4) * nEdgeNodes(9) > 0
    elementFaceNodeIndexes = zeros(1, nEdgeNodes(4) * nEdgeNodes(9));
    c = 0;
    for k = 2:nElementNodesZ - 1
        for j = 2:nElementNodesY - 1
            c = c + 1;
            elementFaceNodeIndexes(c) = (k - 1) * nElementNodesX * nElementNodesY + (j - 1) * nElementNodesX + 1;
        end
    end
    Face(1).elementFaceNodeIndexes = elementFaceNodeIndexes;
else
    Face(1).elementFaceNodeIndexes = [];
end

if nEdgeNodes(2) * nEdgeNodes(10) > 0
    elementFaceNodeIndexes = zeros(1, nEdgeNodes(2) * nEdgeNodes(10));
    c = 0;
    for k = 2:nElementNodesZ - 1
        for j = 2:nElementNodesY - 1
            c = c + 1;
            elementFaceNodeIndexes(c) = (k - 1) * nElementNodesX * nElementNodesY + j * nElementNodesX;
        end
    end
    Face(2).elementFaceNodeIndexes = elementFaceNodeIndexes;
else
    Face(2).elementFaceNodeIndexes = [];
end

if nEdgeNodes(1) * nEdgeNodes(10) > 0
    elementFaceNodeIndexes = zeros(1, nEdgeNodes(1) * nEdgeNodes(10));
    c = 0;
    for k = 2:nElementNodesZ - 1
        for j = 2:nElementNodesX - 1
            c = c + 1;
            elementFaceNodeIndexes(c) = (k - 1) * nElementNodesX * nElementNodesY + j;
        end
    end
    Face(3).elementFaceNodeIndexes = elementFaceNodeIndexes;
else
    Face(3).elementFaceNodeIndexes = [];
end

if nEdgeNodes(3) * nEdgeNodes(11) > 0
    elementFaceNodeIndexes = zeros(1, nEdgeNodes(3) * nEdgeNodes(11));
    c = 0;
    for k = 2:nElementNodesZ - 1
        for j = 2:nElementNodesX - 1
            c = c + 1;
            elementFaceNodeIndexes(c) = (k) * nElementNodesX * nElementNodesY - nElementNodesY + j;
        end
    end
    Face(4).elementFaceNodeIndexes = elementFaceNodeIndexes;
else
    Face(4).elementFaceNodeIndexes = [];
end

if nEdgeNodes(1) * nEdgeNodes(2) > 0
    elementFaceNodeIndexes = zeros(1, nEdgeNodes(1) * nEdgeNodes(2));
    c = 0;
    for k = 2:nElementNodesY - 1
        for j = 2:nElementNodesX - 1
            c = c + 1;
            elementFaceNodeIndexes(c) = (k - 1) * nElementNodesX + j;
        end
    end
    Face(5).elementFaceNodeIndexes = elementFaceNodeIndexes;
else
    Face(5).elementFaceNodeIndexes = [];
end

if nEdgeNodes(5) * nEdgeNodes(6) > 0
    elementFaceNodeIndexes = zeros(1, nEdgeNodes(5) * nEdgeNodes(6));
    c = 0;
    for k = 2:nElementNodesY - 1
        for j = 2:nElementNodesX - 1
            c = c + 1;
            elementFaceNodeIndexes(c) = (k - 1) * nElementNodesX + ...
                nElementNodesX * nElementNodesY * nElementNodesZ - nElementNodesX * nElementNodesY + j;
        end
    end
    Face(6).elementFaceNodeIndexes = elementFaceNodeIndexes;
else
    Face(6).elementFaceNodeIndexes = [];
end

%% interior nodes
if nEdgeNodes(1) * nEdgeNodes(2) * nEdgeNodes(9) > 0
    elementInteriorNodeIndexes = zeros(1, nEdgeNodes(1) * nEdgeNodes(2) * nEdgeNodes(9));
    c = 0;
    for k = 2:nElementNodesZ - 1
        for j = 2:nElementNodesY - 1
            for i = 2:nElementNodesX - 1
                c = c + 1;
                elementInteriorNodeIndexes(c) = (k - 1) * nElementNodesX * nElementNodesY + ...
                                                (j - 1) * nElementNodesY + i;
            end
        end
    end
else
    elementInteriorNodeIndexes = [];
end

% ---------------------------------------------------------------------------------------------------

end
