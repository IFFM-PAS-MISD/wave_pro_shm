function [SemMesh] = gmsh2sem_mesh(msh, PhysicalNames, elementOrder, in_millimeters)
%
% Translate GMSH mesh to Spectral Element Method mesh
% It should be used in conjunction with ``read_mesh.m``
% It reorders nodes in elements
% It redistributes node coordinates according to Gauss-Lobatto-Legendre points
%
% USAGE::
%
%   [SemMesh] = gmsh2sem_mesh(msh,PhysicalNames,elementOrder,in_millimeters)
%
% Arguments:
%     msh (struct):
%       mesh outputed by ``read_msh.m``, structure with fields
%
%       nbNod(integer):
%         number of nodes
%
%       POS (double):
%         coordinates [x,y,z]
%
%       MAX (double):
%         [max(x), max(y), max(z)] coordinates
%
%       MIN (double):
%         [min(x), min(y), min(z)] coordinates
%
%       PNT (double):
%         points, columns (pointIndex, pointTag), dimensions [nPoints,2]
%
%       LINESi (double):
%         lines,i=3,4,5,6, columns (lineNodesList, curveTag)
%         dimensions are order dependent [nLines,(3,4,5,6,7)]
%
%       QUADSi (double):
%         quad elements, i=9,16,25,36,49, columns (element2DList, surfaceTag)
%         dimensions are order dependent [nElements2D, (5,10,17,26,37,50)]
%
%       HEXASi (double):
%         hexahedral elements, i=27,64,125,216,343,
%         columns: (element3DList, volumeTag)
%         dimensions are order dependent [nElements3D, (9,28,65,126,217,344)]
%
%    PhysicalNames (struct):
%      structure outputed by ``read_mesh.m`` containing
%      physical names and tags for parts of the mesh
%
%      dimension (integer):
%        dimension defining points (0), curves (1), surfaces (2) and volumes (3)
%
%      physicalTag (integer):
%        tag, the same as in the last column of elements
%
%      name (string):
%        name describing a part of the mesh
%
%
%    elementOrder (integer):
%      element order in range 1-6
%
%    in_millimeters (logical):
%      flag of model units, true - if mesh prepared in [mm], false - if mesh prepared in [m]
%
% Returns:
%     SemMesh (struct):
%       structure with fields
%
%       elementOrder (integer):
%         element order
%
%       element1DNodes (integer):
%         1D spectral elements topology (element nodes connectivity),
%         dimensions [nElement1D in mesh, nNodes in element1D]
%
%       element2DNodes (integer):
%         2D spectral elements topology (element nodes connectivity),
%         dimensions [nElement2D in mesh, nNodes in element2D]
%
%       element3DNodes (integer):
%         3D spectral elements topology (element nodes connectivity),
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
%         e.g. for material assignment
%
% .. seealso:: Functions :func:`locate_elements_in_physical_entities`,
%                        :func:`gll`,
%                        :func:`nodal_base_change_3d`
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

% extract elements based on PhysicalNames
[SemMesh] = locate_elements_in_physical_entities(msh, PhysicalNames, elementOrder);

% translate node ordering to spectral nodes
% order 1
gmsh2spectralElement1DRenumbering1 = [1, 2];
gmsh2spectralElement2DRenumbering1 = [1, 2, 4, 3];
gmsh2spectralElement3DRenumbering1 = [1, 2, 4, 3, 5, 6, 8, 7];
% order 2
gmsh2spectralElement1DRenumbering2 = [1, 3, 2];
gmsh2spectralElement2DRenumbering2 = [1, 5, 2, 8, 9, 6, 4, 7, 3];
gmsh2spectralElement3DRenumbering2 = [1, 9, 2, 10, 21, 12, 4, 14, 3, 11, 22, 13, 23, 27, 24, 16, 25, 15, 5, 17, 6, 18, 26, 19, 8, 20, 7];
% order 3
gmsh2spectralElement1DRenumbering3 = [1, 3, 4, 2];
gmsh2spectralElement2DRenumbering3 = [1, 5, 6, 2, 12, 13, 14, 7, 11, 16, 15, 8, 4, 10, 9, 3];
gmsh2spectralElement3DRenumbering3 = [1, 9, 10, 2, 11, 33, 36, 15, 12, 34, 35, 16, 4, 20, 19, 3, 13, 37, 38, 17, 41, 57, 58, 45, 44, 60, 59, 46, 23, 50, 49, 21, 14, 40, 39, 18, 42, 61, 62, 48, 43, 64, 63, 47, 24, 51, 52, 22, 5, 25, 26, 6, 27, 53, 54, 29, 28, 56, 55, 30, 8, 32, 31, 7];
% order 4
gmsh2spectralElement1DRenumbering4 = [1, 3, 4, 5, 2];
gmsh2spectralElement2DRenumbering4 = [1, 5, 6, 7, 2, 16, 17, 21, 18, 8, 15, 24, 25, 22, 9, 14, 20, 23, 19, 10, 4, 13, 12, 11, 3];

gmsh2spectralElement3DRenumbering4 = [1, 9, 10, 11, 2, 12, 45, 52, 48, 18, 13, 49, 53, 51, 19, 14, 46, 50, 47, 20, 4, 26, 25, 24, 3, ...
                                      15, 54, 58, 55, 21, 63, 99, 107, 100, 72, 70, 108, 119, 110, 76, 66, 102, 112, 101, 73, 30, 82, 85, 81, 27, ...
                                      16, 61, 62, 59, 22, 67, 109, 120, 111, 79, 71, 121, 125, 122, 80, 69, 114, 123, 113, 77, 31, 86, 89, 88, 28, ...
                                      17, 57, 60, 56, 23, 64, 103, 115, 104, 75, 68, 116, 124, 117, 78, 65, 106, 118, 105, 74, 32, 83, 87, 84, 29, ...
                                      5, 33, 34, 35, 6, 36, 90, 94, 91, 39, 37, 97, 98, 95, 40, 38, 93, 96, 92, 41, 8, 44, 43, 42, 7];
% order 5
gmsh2spectralElement1DRenumbering5 = [1, 3, 4, 5, 6, 2];
gmsh2spectralElement2DRenumbering5 = [1, 5, 6, 7, 8, 2, 20, 21, 25, 22, 9, 19, 32, 33, 34, 27, 10, 18, 31, 36, 35, 28, 11, 17, 24, 30, 29, 23, 12, 4, 16, 15, 14, 13, 3];

gmsh2spectralElement3DRenumbering5 = [1, 9:12, 2, 13, 57, 68, 67, 60, 21, 14, 61, 69, 72, 66, 22, 15, 62, 70, 71, 65, 23, 16, 58, 63, 64, 59, 24, 4, 32, 31, 30, 29, 3, ...
                                      17, 73, 77, 78, 74, 25, 89, 153, 161, 162, 154, 105, 100, 163, 185, 188, 167, 109, 99, 164, 186, 187, 168, 110, 92, 156, 172, 171, 155, 106, 37, 122, 126, 125, 121, 33, ...
                                      18, 84, 85, 86, 79, 26, 93, 165, 189, 190, 169, 116, 101, 193, 209, 210, 197, 117, 104, 196, 212, 211, 198, 118, 98, 175, 202, 201, 173, 111, 38, 127, 134, 133, 132, 34, ...
                                      19, 83, 88, 87, 80, 27, 94, 166, 192, 191, 170, 115, 102, 194, 213, 214, 200, 120, 103, 195, 216, 215, 199, 119, 97, 176, 203, 204, 174, 112, 39, 128, 135, 136, 131, 35, ...
                                      20, 76, 82, 81, 75, 28, 90, 157, 177, 178, 158, 108, 95, 179, 205, 206, 181, 114, 96, 180, 208, 207, 182, 113, 91, 160, 184, 183, 159, 107, 40, 123, 129, 130, 124, 36, ...
                                      5, 41, 42, 43, 44, 6, 45, 137, 141, 142, 138, 49, 46, 148, 149, 150, 143, 50, 47, 147, 152, 151, 144, 51, 48, 140, 146, 145, 139, 52, 8, 56, 55, 54, 53, 7];

% order 6
gmsh2spectralElement1DRenumbering6 = [1, 3, 4, 5, 6, 7, 2];
gmsh2spectralElement2DRenumbering6 = [1, 5:9, 2, 24, 25, 29, 30, 31, 26, 10, 23, 40, 41, 45, 42, 32, 11, 22, 39, 48, 49, 46, 33, 12, 21, 38, 44, 47, 43, 34, 13, 20, 28, 37, 36, 35, 27, 14, 4, 19, 18, 17, 16, 15, 3];

gmsh2spectralElement3DRenumbering6 = [1, 9:13, 2, 14, 69, 84, 83, 82, 72, 24, 15, 73, 85, 92, 88, 81, 25, 16, 74, 89, 93, 91, 80, 26, 17, 75, 86, 90, 87, 79, 27, 18, 70, 76, 77, 78, 71, 28, 4, 38, 37, 36, 35, 34, 3, ...
                                      19, 94, 98, 99, 100, 95, 29, 119, 219, 227, 228, 229, 220, 144, 134, 230, 263, 270, 266, 236, 148, 133, 231, 267, 271, 269, 237, 149, 132, 232, 264, 268, 265, 238, 150, 122, 222, 244, 243, 242, 221, 145, 44, 170, 175, 174, 173, 169, 39, ...
                                      20, 109, 110, 114, 111, 101, 30, 123, 233, 272, 276, 273, 239, 159, 135, 281, 317, 325, 318, 290, 160, 142, 288, 326, 337, 328, 294, 164, 138, 284, 320, 330, 319, 291, 161, 131, 248, 300, 303, 299, 245, 151, 45, 176, 186, 189, 185, 184, 40, ...
                                      21, 108, 117, 118, 115, 102, 31, 124, 234, 279, 280, 277, 240, 158, 139, 285, 327, 338, 329, 297, 167, 143, 289, 339, 343, 340, 298, 168, 141, 287, 332, 341, 331, 295, 165, 130, 249, 304, 307, 306, 246, 152, 46, 177, 190, 193, 192, 183, 41, ...
                                      22, 107, 113, 116, 112, 103, 32, 125, 235, 275, 278, 274, 241, 157, 136, 282, 321, 333, 322, 293, 163, 140, 286, 334, 342, 335, 296, 166, 137, 283, 324, 336, 323, 292, 162, 129, 250, 301, 305, 302, 247, 153, 47, 178, 187, 191, 188, 182, 42, ...
                                      23, 97, 106, 105, 104, 96, 33, 120, 223, 251, 252, 253, 224, 147, 126, 254, 308, 312, 309, 257, 156, 127, 255, 315, 316, 313, 258, 155, 128, 256, 311, 314, 310, 259, 154, 121, 226, 262, 261, 260, 225, 146, 48, 171, 179, 180, 181, 172, 43, ...
                                      5, 49, 50, 51, 52, 53, 6, 54, 194, 198, 199, 200, 195, 59, 55, 209, 210, 214, 211, 201, 60, 56, 208, 217, 218, 215, 202, 61, 57, 207, 213, 216, 212, 203, 62, 58, 197, 206, 205, 204, 196, 63, 8, 68, 67, 66, 65, 64, 7];

if elementOrder == 1
    if isfield(msh, 'LINES')
        SemMesh.element1DNodes = msh.LINES(:, gmsh2spectralElement1DRenumbering1);
    end
    if isfield(msh, 'QUADS')
        SemMesh.element2DNodes = msh.QUADS(:, gmsh2spectralElement2DRenumbering1);
    end
    if isfield(msh, 'HEXAS')
        SemMesh.element3DNodes = msh.HEXAS(:, gmsh2spectralElement3DRenumbering1);
    end
end

if elementOrder == 2
    if isfield(msh, 'LINES3')
        SemMesh.element1DNodes = msh.LINES3(:, gmsh2spectralElement1DRenumbering2);
    end
    if isfield(msh, 'QUADS9')
        SemMesh.element2DNodes = msh.QUADS9(:, gmsh2spectralElement2DRenumbering2);
    end
    if isfield(msh, 'HEXAS27')
        SemMesh.element3DNodes = msh.HEXAS27(:, gmsh2spectralElement3DRenumbering2);
    end
end

if elementOrder == 3
    if isfield(msh, 'LINES4')
        SemMesh.element1DNodes = msh.LINES4(:, gmsh2spectralElement1DRenumbering3);
    end
    if isfield(msh, 'QUADS16')
        SemMesh.element2DNodes = msh.QUADS16(:, gmsh2spectralElement2DRenumbering3);
    end
    if isfield(msh, 'HEXAS64')
        SemMesh.element3DNodes = msh.HEXAS64(:, gmsh2spectralElement3DRenumbering3);
    end
end

if elementOrder == 4
    if isfield(msh, 'LINES5')
        SemMesh.element1DNodes = msh.LINES5(:, gmsh2spectralElement1DRenumbering4);
    end
    if isfield(msh, 'QUADS25')
        SemMesh.element2DNodes = msh.QUADS25(:, gmsh2spectralElement2DRenumbering4);
    end
    if isfield(msh, 'HEXAS125')
        SemMesh.element3DNodes = msh.HEXAS125(:, gmsh2spectralElement3DRenumbering4);
    end
end

if elementOrder == 5
    if isfield(msh, 'LINES6')
        SemMesh.element1DNodes = msh.LINES6(:, gmsh2spectralElement1DRenumbering5);
    end
    if isfield(msh, 'QUADS36')
        SemMesh.element2DNodes = msh.QUADS36(:, gmsh2spectralElement2DRenumbering5);
    end
    if isfield(msh, 'HEXAS216')
        SemMesh.element3DNodes = msh.HEXAS216(:, gmsh2spectralElement3DRenumbering5);
    end
end

if elementOrder == 6
    if isfield(msh, 'LINES7')
        SemMesh.element1DNodes = msh.LINES7(:, gmsh2spectralElement1DRenumbering6);
    end
    if isfield(msh, 'QUADS49')
        SemMesh.element2DNodes = msh.QUADS49(:, gmsh2spectralElement2DRenumbering6);
    end
    if isfield(msh, 'HEXAS349')
        SemMesh.element3DNodes = msh.HEXAS343(:, gmsh2spectralElement3DRenumbering6);
    end
end

% scale coordinates if needed to units in [m]
if in_millimeters
    SemMesh.nodeCoordinates = msh.POS / 1e3;
else
    SemMesh.nodeCoordinates = msh.POS;
end

% adjust coordinates according to Gauss-Lobatto-Legendre points
[ksi, ~] = gll(elementOrder + 1); % Gauss-Lobatto-Legendre points
eta = ksi;
dzeta = ksi;
ksiOld = linspace(-1, 1, elementOrder + 1);
etaOld = linspace(-1, 1, elementOrder + 1);
dzetaOld = linspace(-1, 1, elementOrder + 1);
newCoords = zeros(size(SemMesh.nodeCoordinates));
nElments3D = size(SemMesh.element3DNodes, 1);
for iElement = 1:nElments3D
    % high order interpolation of geometry
    elementNodes = SemMesh.element3DNodes(iElement, :);
    elementCoordsX = SemMesh.nodeCoordinates(elementNodes, 1)';
    elementCoordsY = SemMesh.nodeCoordinates(elementNodes, 2)';
    elementCoordsZ = SemMesh.nodeCoordinates(elementNodes, 3)';
    [xNew, yNew, zNew] = nodal_base_change_3d(ksi, eta, dzeta, ksiOld, etaOld, dzetaOld, elementCoordsX, elementCoordsY, elementCoordsZ);

    newCoords(elementNodes, 1) = xNew;
    newCoords(elementNodes, 2) = yNew;
    newCoords(elementNodes, 3) = zNew;
end
SemMesh.nodeCoordinates = newCoords;
% ---------------------------------------------------------------------------------------------------

end
