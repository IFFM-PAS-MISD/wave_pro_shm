function [SemMesh] = locate_elements_in_physical_entities(msh, PhysicalNames, elementOrder)
%
% Locate elements in physical entities.
% Regroup physical entities based on dimensionality and Physical Names
%
% USAGE::
%
%   [SemMesh] = locate_elements_in_physical_entities(msh,PhysicalNames,elementOrder)
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
%         [min(x), min(y) min(z)] coordinates
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
%      structure outputed by read_mesh.m containing
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
%    elementOrder (integer):
%      element order in range 1-6
%
% Returns:
%     SemMesh (struct):
%       structure with fields
%
%       elementOrder (integer):
%         element order
%
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
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

SemMesh.elementOrder = elementOrder;

% extract elements based on PhysicalNames
PhysicalVolumes_idx = 0;
PhysicalPoints_idx = 0;
PhysicalSurfaces_idx = 0;
PhysicalLines_idx = 0;
for i = 1:length(PhysicalNames)
    if PhysicalNames(i).dimension == 0
        PhysicalPoints_idx = PhysicalPoints_idx + 1;
        SemMesh.PhysicalPoints(PhysicalPoints_idx).nodes = ...
            msh.PNT((find(msh.PNT(:, 2) == PhysicalNames(i).physicalTag))', 1);
        SemMesh.PhysicalPoints(PhysicalPoints_idx).physicalTag = PhysicalNames(i).physicalTag;
        SemMesh.PhysicalPoints(PhysicalPoints_idx).name = PhysicalNames(i).name;
    end
    if PhysicalNames(i).dimension == 1
        PhysicalLines_idx = PhysicalLines_idx + 1;
        SemMesh.PhysicalLines(PhysicalLines_idx).physicalTag = PhysicalNames(i).physicalTag;
        SemMesh.PhysicalLines(PhysicalLines_idx).name = PhysicalNames(i).name;
        if elementOrder == 1
            indexes = (find(msh.LINES(:, end) == PhysicalNames(i).physicalTag))';
            SemMesh.PhysicalLines(PhysicalLines_idx).nodes = msh.LINES(indexes, 1:end - 1);
        end
        if elementOrder == 2
            indexes = (find(msh.LINES3(:, end) == PhysicalNames(i).physicalTag))';
            SemMesh.PhysicalLines(PhysicalLines_idx).nodes = msh.LINES3(indexes, 1:end - 1);
        end
        if elementOrder == 3
            indexes = (find(msh.LINES4(:, end) == PhysicalNames(i).physicalTag))';
            SemMesh.PhysicalLines(PhysicalLines_idx).nodes = msh.LINES4(indexes, 1:end - 1);
        end
        if elementOrder == 4
            indexes = (find(msh.LINES5(:, end) == PhysicalNames(i).physicalTag))';
            SemMesh.PhysicalLines(PhysicalLines_idx).nodes = msh.LINES5(indexes, 1:end - 1);
        end
        if elementOrder == 5
            indexes = (find(msh.LINES6(:, end) == PhysicalNames(i).physicalTag))';
            SemMesh.PhysicalLines(PhysicalLines_idx).nodes = msh.LINES6(indexes, 1:end - 1);
        end
        if elementOrder == 6
            indexes = (find(msh.LINES7(:, end) == PhysicalNames(i).physicalTag))';
            SemMesh.PhysicalLines(PhysicalLines_idx).nodes = msh.LINES7(indexes, 1:end - 1);
        end
    end
    if PhysicalNames(i).dimension == 2
        PhysicalSurfaces_idx = PhysicalSurfaces_idx + 1;
        SemMesh.PhysicalSurfaces(PhysicalSurfaces_idx).physicalTag = PhysicalNames(i).physicalTag;
        SemMesh.PhysicalSurfaces(PhysicalSurfaces_idx).name = PhysicalNames(i).name;
        if elementOrder == 1
            SemMesh.PhysicalSurfaces(PhysicalSurfaces_idx).elements = ...
                (find(msh.QUADS(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 2
            SemMesh.PhysicalSurfaces(PhysicalSurfaces_idx).elements = ...
                (find(msh.QUADS9(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 3
            SemMesh.PhysicalSurfaces(PhysicalSurfaces_idx).elements = ...
                (find(msh.QUADS16(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 4
            SemMesh.PhysicalSurfaces(PhysicalSurfaces_idx).elements = ...
                (find(msh.QUADS25(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 5
            SemMesh.PhysicalSurfaces(PhysicalSurfaces_idx).elements = ...
                (find(msh.QUADS36(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 6
            SemMesh.PhysicalSurfaces(PhysicalSurfaces_idx).elements = ...
                (find(msh.QUADS49(:, end) == PhysicalNames(i).physicalTag))';
        end
    end
    if PhysicalNames(i).dimension == 3
        PhysicalVolumes_idx = PhysicalVolumes_idx + 1;
        SemMesh.PhysicalVolumes(PhysicalVolumes_idx).physicalTag = PhysicalNames(i).physicalTag;
        SemMesh.PhysicalVolumes(PhysicalVolumes_idx).name = PhysicalNames(i).name;
        if elementOrder == 1
            SemMesh.PhysicalVolumes(PhysicalVolumes_idx).elements = ...
                (find(msh.HEXAS(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 2
            SemMesh.PhysicalVolumes(PhysicalVolumes_idx).elements = ...
                (find(msh.HEXAS27(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 3
            SemMesh.PhysicalVolumes(PhysicalVolumes_idx).elements = ...
                (find(msh.HEXAS64(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 4
            SemMesh.PhysicalVolumes(PhysicalVolumes_idx).elements = ...
                (find(msh.HEXAS125(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 5
            SemMesh.PhysicalVolumes(PhysicalVolumes_idx).elements = ...
                (find(msh.HEXAS216(:, end) == PhysicalNames(i).physicalTag))';
        end
        if elementOrder == 6
            SemMesh.PhysicalVolumes(PhysicalVolumes_idx).elements = ...
                (find(msh.HEXAS343(:, end) == PhysicalNames(i).physicalTag))';
        end
    end
end

% ---------------------------------------------------------------------------------------------------

end
