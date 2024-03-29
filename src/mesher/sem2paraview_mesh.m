function [ParaviewMesh] = sem2paraview_mesh(SemMesh)
%
% Translate Spectral Element Method mesh to Paraview (modifies connectivity in element nodes)
%
% USAGE::
%
%   [ParaviewMesh] = sem2paraview_mesh(SemMesh)
%
% Arguments:
%     SemMesh (struct):
%       Spectral Element Method mesh with fields
%
%       elementOrder (integer):
%         element order
%
%       element1DNodes (integer):
%         1D spectral elements topology (element nodes connectivity)
%         dimensions [nElement1D in mesh, nNodes in element1D]
%
%       element2DNodes (integer):
%         2D spectral elements topology (element nodes connectivity)
%         dimensions [nElement2D in mesh, nNodes in element2D]
%
%       element3DNodes (integer):
%         3D spectral elements topology (element nodes connectivity)
%         dimensions [nElement3D in mesh, nNodes in element3D]
%
%       nodeCoordinates (double):
%         coordinates of all mesh nodes [x,y,z],
%         dimensions [nNodes,3], Units [m]
%
%       physicalPoints (struct):
%         set of points and corresponding physical names,
%         e.g. for point force or outputs
%
%       physicalLines (struct):
%         set of lines and corresponding physical names,
%         e.g. for damage boundary
%
%       physicalSurfaces (struct):
%         set of surfaces and corresponding physical names,
%         e.g. for boundary conditions, piezoelectric electrodes
%
%       physicalVolumes (struct):
%         set of volumes and corresponding physical names,
%         e.g. for material assigment
%
% Returns:
%     ParaviewMesh (struct):
%       Parview mesh consisting of LAGRANGEHEXAHEDRON elements with fields
%
%       element1DNodes (integer):
%         1D spectral elements topology (element nodes connectivity)
%         dimensions [nElement1D in mesh, nNodes in element1D]
%
%       element2DNodes (integer):
%         2D spectral elements topology (element nodes connectivity)
%         dimensions [nElement2D in mesh, nNodes in element2D]
%
%       element3DNodes (integer):
%         3D spectral elements topology (element nodes connectivity)
%         dimensions [nElement3D in mesh, nNodes in element3D]
%
%       nodeCoordinates (double):
%         coordinates of all mesh nodes [x,y,z],
%         dimensions [nNodes,3], Units [m]
%
%
% .. seealso:: Function :func:`element3d_numbering`
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

elementOrder = SemMesh.elementOrder;
nElementNodesX = elementOrder + 1;
nElementNodesY = elementOrder + 1;
nElementNodesZ = elementOrder + 1;

[elementCornerNodeIndexes, Edge, Face, elementInteriorNodeIndexes] = ...
 element3d_numbering(nElementNodesX, nElementNodesY, nElementNodesZ);
 
elementEdgeNodeIndexes = [Edge(1).elementEdgeNodeIndexes, Edge(2).elementEdgeNodeIndexes, ...
                          Edge(3).elementEdgeNodeIndexes, Edge(4).elementEdgeNodeIndexes, ...
                          Edge(5).elementEdgeNodeIndexes, Edge(6).elementEdgeNodeIndexes, ...
                          Edge(7).elementEdgeNodeIndexes, Edge(8).elementEdgeNodeIndexes, ...
                          Edge(9).elementEdgeNodeIndexes, Edge(10).elementEdgeNodeIndexes, ...
                          Edge(11).elementEdgeNodeIndexes, Edge(12).elementEdgeNodeIndexes];
elementFaceNodeIndexes = [Face(1).elementFaceNodeIndexes, Face(2).elementFaceNodeIndexes, ...
                          Face(3).elementFaceNodeIndexes, Face(4).elementFaceNodeIndexes, ...
                          Face(5).elementFaceNodeIndexes, Face(6).elementFaceNodeIndexes];
% translate node ordering to spectral nodes
% order 1
spectralElement1D2ParaviewRenumbering1 = [1, 2];
spectralElement2D2ParaviewRenumbering1 = [1, 2, 4, 3];
spectralElement3D2ParaviewRenumbering1 = [1, 2, 4, 3, 5, 6, 8, 7];
% order 2
spectralElement1D2ParaviewRenumbering2 = [1, 3, 2];
spectralElement2D2ParaviewRenumbering2 = [1, 3, 9, 7, 2, 6, 8, 4, 5];
spectralElement3D2ParaviewRenumbering2 = [1, 3, 9, 7, 19, 21, 27, 25, 2, 6, 8, 4, 20, 24, 26, 22, 10, 12, 18, 16, 13, 15, 11, 17, 5, 23, 14];
% order 3-6
spectralElement1D2ParaviewRenumbering = [1, elementCornerNodeIndexes(2), Edge(1).elementEdgeNodeIndexes];
spectralElement2D2ParaviewRenumbering = [elementCornerNodeIndexes([1, 2, 4, 3]), ...
                                         Edge(1).elementEdgeNodeIndexes, ...
                                         Edge(2).elementEdgeNodeIndexes, ...
                                         Edge(3).elementEdgeNodeIndexes, ...
                                         Edge(4).elementEdgeNodeIndexes, ...
                                         Face(5).elementFaceNodeIndexes];
spectralElement3D2ParaviewRenumbering = [elementCornerNodeIndexes([1, 2, 4, 3, 5, 6, 8, 7]), ...
                                         elementEdgeNodeIndexes, elementFaceNodeIndexes, elementInteriorNodeIndexes];

if elementOrder == 1
    if isfield(SemMesh, 'element1DNodes')
        ParaviewMesh.element1DNodes = SemMesh.element1DNodes(:, spectralElement1D2ParaviewRenumbering1) - 1;
    end
    if isfield(SemMesh, 'element2DNodes')
        ParaviewMesh.element2DNodes = SemMesh.element2DNodes(:, spectralElement2D2ParaviewRenumbering1) - 1;
    end
    if isfield(SemMesh, 'element3DNodes')
        ParaviewMesh.element3DNodes = SemMesh.element3DNodes(:, spectralElement3D2ParaviewRenumbering1) - 1;
    end
end

if elementOrder == 2
    if isfield(SemMesh, 'element1DNodes')
        ParaviewMesh.element1DNodes = SemMesh.element1DNodes(:, spectralElement1D2ParaviewRenumbering2) - 1;
    end
    if isfield(SemMesh, 'element2DNodes')
        ParaviewMesh.element2DNodes = SemMesh.element2DNodes(:, spectralElement2D2ParaviewRenumbering2) - 1;
    end
    if isfield(SemMesh, 'element3DNodes')
        ParaviewMesh.element3DNodes = SemMesh.element3DNodes(:, spectralElement3D2ParaviewRenumbering2) - 1;
    end
end

if elementOrder >= 3
    if isfield(SemMesh, 'element1DNodes')
        ParaviewMesh.element1DNodes = SemMesh.element1DNodes(:, spectralElement1D2ParaviewRenumbering) - 1;
    end
    if isfield(SemMesh, 'element2DNodes')
        ParaviewMesh.element2DNodes = SemMesh.element2DNodes(:, spectralElement2D2ParaviewRenumbering) - 1;
    end
    if isfield(SemMesh, 'element3DNodes')
        ParaviewMesh.element3DNodes = SemMesh.element3DNodes(:, spectralElement3D2ParaviewRenumbering) - 1;
    end
end

ParaviewMesh.nodeCoordinates = SemMesh.nodeCoordinates;

% ---------------------------------------------------------------------------------------------------

end
