function [msh, PhysicalNames] = read_msh_mesh(meshfile, isVerbose)
%
% Read .msh 4.1 file format from `GMSH <https://gmsh.info/>`_ and save as .mat file
% The .msh file format description: http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
% read_msh_mesh is focused on reading hexahedral complete elements of orders 1-6
% and corresponding nodes and PhysicalNames if they exist.
% It was created because large .m files with meshes exported from GMSH
% cannot be loaded to MATLAB by using 'run' command.
% See also python gmshparser for alternative implementation:
% https://gmshparser.readthedocs.io/
%
% USAGE::
%
%   [msh,PhysicalNames] = read_msh_mesh(meshfile,isVerbose)
%
% Arguments:
%     meshfile (string):
%       filename including path and extension
%
%     isVerbose (logical):
%       flag for spitting additional information during calculations
%
% Returns:
%     msh (struct):
%       mesh outputed by read_msh.m, structure with fields
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
%
%    elementOrder (integer):
%      element order in range 1-6
%
%    in_millimeters (logical):
%      flag of model units, true - if mesh prepared in [mm], false - if mesh prepared in [m]
%
%
% .. Note:: Variable names are kept the same as in the description of mesh
%       msh format version 4 (current revision: version 4.1)
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

if nargin < 2
    output_path = '';
end

EntitiesExist = false;
PhysicalNamesExist = false;

[filepath, mesh_name, ext] = fileparts(meshfile);

% check total number of lines
fid = fopen(meshfile);
if fid == -1
    if isVerbose
        disp(['error: file ', meshfile, ' does not exist!']);
    end
    return
end
line_counter = 0;
while 1
    line_counter = line_counter + 1;
    str = fgetl(fid);
    if str == -1
        total_lines = line_counter;
        break
    end
end

% $PhysicalNames
frewind(fid);
line_counter = 0;
while 1
    line_counter = line_counter + 1;
    str = fgetl(fid);
    TF = startsWith(str, '$PhysicalNames');
    if TF == true
        break
    end
end

% $PhysicalNames are optional so check if they exist
if line_counter < total_lines
    PhysicalNamesExist = true;
end
if PhysicalNamesExist
    if isVerbose
        disp('reading PhysicalNames');
    end
    str = fgetl(fid);
    C = textscan(str, '%d', 1, 'delimiter', '\n');
    line_counter = line_counter + 1;
    numPhysicalNames = C{1};
    PhysicalNames(numPhysicalNames) = struct();
    for i = 1:numPhysicalNames
        str = fgetl(fid);
        C = textscan(str, '%d %d %s', 1, 'delimiter', '\n');
        % C = textscan(fid,'%d %d %s',1,'delimiter','\n', 'headerlines',line_counter); % alternative way
        PhysicalNames(i).dimension = C{1};
        PhysicalNames(i).physicalTag = C{2};
        temp = char(C{3});
        PhysicalNames(i).name = temp(2:end - 1);

        line_counter = line_counter + 1;
    end
end
% $Entities
frewind(fid);
line_counter = 0;
while 1
    line_counter = line_counter + 1;
    str = fgetl(fid);
    TF = startsWith(str, '$Entities');
    if TF == true
        break
    end
end
% $Entities are optional so check if they exist
if line_counter < total_lines
    EntitiesExist = true;
end
if EntitiesExist
    if isVerbose
        disp('reading Entities');
    end
    str = fgetl(fid);
    C = textscan(str, '%d %d %d %d', 1, 'delimiter', '\n');
    numPoints = C{1};
    numCurves = C{2};
    numSurfaces = C{3};
    numVolumes = C{4};
    for i = 1:numPoints
        str = fgetl(fid);
        C = textscan(str, '%d %d %d %d %d %d', 1, 'delimiter', '\n');
        pointTag = C{1};
        if C{5} > 0  % PhysicalTag exist
            physicalTag.point{pointTag} = C{6};
        end
    end
    for i = 1:numCurves
        str = fgetl(fid);
        C = textscan(str, '%d %d %d %d %d %d %d %d %d', 1, 'delimiter', '\n');
        curveTag = C{1};
        if C{8} > 0  % PhysicalTag exist
            physicalTag.curve{curveTag} = C{9};
        end
    end
    for i = 1:numSurfaces
        str = fgetl(fid);
        C = textscan(str, '%d %d %d %d %d %d %d %d %d', 1, 'delimiter', '\n');
        surfaceTag = C{1};
        if C{8} > 0  % PhysicalTag exist
            physicalTag.surface{surfaceTag} = C{9};
        end
    end
    for i = 1:numVolumes
        str = fgetl(fid);
        C = textscan(str, '%d %d %d %d %d %d %d %d %d', 1, 'delimiter', '\n');
        volumeTag = C{1};
        if C{8} > 0  % PhysicalTag exist
            physicalTag.volume{volumeTag} = C{9};
        end
    end
end

% $PartitionedEntities - are skipped

% $Nodes
if isVerbose
    disp('reading Nodes');
end
frewind(fid);
line_counter = 0;
while 1
    line_counter = line_counter + 1;
    str = fgetl(fid);
    TF = startsWith(str, '$Nodes');
    if TF == true
        break
    end
end
% nodes_line_no = line_counter; % for debugging
str = fgetl(fid);
line_counter = line_counter + 1;
C = textscan(str, '%d %d %d %d', 1, 'delimiter', '\n');
numEntityBlocks = C{1};
numNodes = C{2};
minNodeTag = C{3};
maxNodeTag = C{4};

coords = zeros(numNodes, 3);

coord_num = 0;
for i = 1:numEntityBlocks
    str = fgetl(fid);
    line_counter = line_counter + 1;
    C = textscan(str, '%d %d %d %d', 1, 'delimiter', '\n');
    entityDim = C{1};
    entityTag = C{2};
    parametric = C{3};
    numNodesInBlock = C{4};
    for j = 1:numNodesInBlock
        str = fgetl(fid);
        % node tags are not saved so we skip these lines - use 2 lines below if you need them
        % C = textscan(str,'%d',1,'delimiter','\n');
        % nodeTag = C{1};
        line_counter = line_counter + 1;
    end

    for j = 1:numNodesInBlock
        str = fgetl(fid);
        line_counter = line_counter + 1;
        coord_num = coord_num + 1;
        C = textscan(str, '%f %f %f', 1, 'delimiter', '\n');
        coords(coord_num, :) = [C{1}, C{2}, C{3}];
    end
end

% $Elements
if isVerbose
    disp('reading Elements');
end
frewind(fid);
line_counter = 0;
while 1
    line_counter = line_counter + 1;
    str = fgetl(fid);
    TF = startsWith(str, '$Elements');
    if TF == true
        break
    end
end
% elements_line_no = line_counter; % for debugging
str = fgetl(fid);
line_counter = line_counter + 1;
C = textscan(str, '%d %d %d %d', 1, 'delimiter', '\n');
numEntityBlocks = C{1};
numElements = C{2};
minElementTag = C{3};
maxElementTag = C{4};
% initialize indexes
points_idx = 0;
lines_idx = 0;
lines3_idx = 0;
lines4_idx = 0;
lines5_idx = 0;
lines6_idx = 0;
lines7_idx = 0;
quads_idx = 0;
quads9_idx = 0;
quads16_idx = 0;
quads25_idx = 0;
quads36_idx = 0;
quads49_idx = 0;
hexas_idx = 0;
hexas27_idx = 0;
hexas64_idx = 0;
hexas125_idx = 0;
hexas216_idx = 0;
hexas343_idx = 0;
for i = 1:numEntityBlocks
    str = fgetl(fid);
    line_counter = line_counter + 1;
    C = textscan(str, '%d %d %d %d', 1, 'delimiter', '\n');
    entityDim = C{1};
    entityTag = C{2};
    elementType = C{3};
    numElementsInBlock = C{4};
    if elementType == 15  % Points
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            points_idx = points_idx + 1;
            C = textscan(str, '%d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = C{2};
            if PhysicalNamesExist && EntitiesExist
                points{points_idx, :} = [nodeTag, physicalTag.point{entityTag}];
            else
                points{points_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 1  % Lines
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            lines_idx = lines_idx + 1;
            C = textscan(str, '%d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:3));
            if PhysicalNamesExist && EntitiesExist
                lines{lines_idx, :} = [nodeTag, physicalTag.curve{entityTag}];
            else
                lines{lines_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 8  % 3-node second order line
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            lines3_idx = lines3_idx + 1;
            C = textscan(str, '%d %d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:4));
            if PhysicalNamesExist && EntitiesExist
                lines3{lines3_idx, :} = [nodeTag, physicalTag.curve{entityTag}];
            else
                lines3{lines3_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 26  % 4-node third order line
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            lines4_idx = lines4_idx + 1;
            C = textscan(str, '%d %d %d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:5));
            if PhysicalNamesExist && EntitiesExist
                lines4{lines4_idx, :} = [nodeTag, physicalTag.curve{entityTag}];
            else
                lines4{lines4_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 27  % 5-node fourth order line
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            lines5_idx = lines5_idx + 1;
            C = textscan(str, '%d %d %d %d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:6));
            if PhysicalNamesExist && EntitiesExist
                lines5{lines5_idx, :} = [nodeTag, physicalTag.curve{entityTag}];
            else
                lines5{lines5_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 28  % 6-node fifth order line
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            lines6_idx = lines6_idx + 1;
            C = textscan(str, '%d %d %d %d %d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:7));
            if PhysicalNamesExist && EntitiesExist
                lines6{lines6_idx, :} = [nodeTag, physicalTag.curve{entityTag}];
            else
                lines6{lines6_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 62  % 7-node sixth order line
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            lines7_idx = lines7_idx + 1;
            C = textscan(str, '%d %d %d %d %d %d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:8));
            if PhysicalNamesExist && EntitiesExist
                lines7{lines7_idx, :} = [nodeTag, physicalTag.curve{entityTag}];
            else
                lines7{lines7_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 3  % 4-node quadrangle
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            quads_idx = quads_idx + 1;
            C = textscan(str, '%d %d %d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:5));
            if PhysicalNamesExist && EntitiesExist
                quads{quads_idx, :} = [nodeTag, physicalTag.surface{entityTag}];
            else
                quads{quads_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 10  % 9-node second order quadrangle
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            quads9_idx = quads9_idx + 1;
            C = textscan(str, '%d %d %d %d %d %d %d %d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:10));
            if PhysicalNamesExist && EntitiesExist
                quads9{quads9_idx, :} = [nodeTag, physicalTag.surface{entityTag}];
            else
                quads9{quads9_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 36  % 16-node third order quadrangle
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            quads16_idx = quads16_idx + 1;
            C = cell2mat(textscan(str, '%d', 17));
            elementTag = C(1);
            nodeTag = C(2:17)';
            if PhysicalNamesExist && EntitiesExist
                quads16{quads16_idx, :} = [nodeTag, physicalTag.surface{entityTag}];
            else
                quads16{quads16_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 37  % 25-node fourth order quadrangle
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            quads25_idx = quads25_idx + 1;
            C = cell2mat(textscan(str, '%d', 26));
            elementTag = C(1);
            nodeTag = C(2:26)';
            if PhysicalNamesExist && EntitiesExist
                quads25{quads25_idx, :} = [nodeTag, physicalTag.surface{entityTag}];
            else
                quads25{quads25_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 38  % 36-node fifth order quadrangle
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            quads36_idx = quads36_idx + 1;
            C = cell2mat(textscan(str, '%d', 37));
            elementTag = C(1);
            nodeTag = C(2:37)';
            if PhysicalNamesExist && EntitiesExist
                quads36{quads36_idx, :} = [nodeTag, physicalTag.surface{entityTag}];
            else
                quads36{quads36_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 47  % 49-node sixth order quadrangle
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            quads49_idx = quads49_idx + 1;
            C = cell2mat(textscan(str, '%d', 50));
            elementTag = C(1);
            nodeTag = C(2:50)';
            if PhysicalNamesExist && EntitiesExist
                quads49{quads49_idx, :} = [nodeTag, physicalTag.surface{entityTag}];
            else
                quads49{quads49_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 5  % 8-node hexahedron
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            hexas_idx = hexas_idx + 1;
            C = textscan(str, '%d %d %d %d %d %d %d %d %d', 1, 'delimiter', '\n');
            elementTag = C{1};
            nodeTag = cell2mat(C(2:9));
            if PhysicalNamesExist && EntitiesExist
                hexas{hexas_idx, :} = [nodeTag, physicalTag.volume{entityTag}];
            else
                hexas{hexas_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 12  % 27-node second order hexahedron
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            hexas27_idx = hexas27_idx + 1;
            C = cell2mat(textscan(str, '%d', 28));
            elementTag = C(1);
            nodeTag = C(2:28)';
            if PhysicalNamesExist && EntitiesExist
                hexas27{hexas27_idx, :} = [nodeTag, physicalTag.volume{entityTag}];
            else
                hexas27{hexas27_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 92  % 64-node third order hexahedron
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            hexas64_idx = hexas64_idx + 1;
            C = cell2mat(textscan(str, '%d', 65));
            elementTag = C(1);
            nodeTag = C(2:65)';
            if PhysicalNamesExist && EntitiesExist
                hexas64{hexas64_idx, :} = [nodeTag, physicalTag.volume{entityTag}];
            else
                hexas64{hexas64_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 93  % 125-node fourth order hexahedron
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            hexas125_idx = hexas125_idx + 1;
            C = cell2mat(textscan(str, '%d', 126));
            elementTag = C(1);
            nodeTag = C(2:126)';
            if PhysicalNamesExist && EntitiesExist
                hexas125{hexas125_idx, :} = [nodeTag, physicalTag.volume{entityTag}];
            else
                hexas125{hexas125_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 94  % 216-node fifth order hexahedron
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            hexas216_idx = hexas216_idx + 1;
            C = cell2mat(textscan(str, '%d', 217));
            elementTag = C(1);
            nodeTag = C(2:217)';
            if PhysicalNamesExist && EntitiesExist
                hexas216{hexas216_idx, :} = [nodeTag, physicalTag.volume{entityTag}];
            else
                hexas216{hexas216_idx, :} = nodeTag;
            end
        end
    end
    if elementType == 95  % 343-node sixth order hexahedron
        for j = 1:numElementsInBlock
            str = fgetl(fid);
            line_counter = line_counter + 1;
            hexas343_idx = hexas343_idx + 1;
            C = cell2mat(textscan(str, '%d', 344));
            elementTag = C(1);
            nodeTag = C(2:344)';
            if PhysicalNamesExist && EntitiesExist
                hexas343{hexas343_idx, :} = [nodeTag, physicalTag.volume{entityTag}];
            else
                hexas343{hexas343_idx, :} = nodeTag;
            end
        end
    end
end
fclose(fid);

msh.nbNod = numNodes;
msh.POS = coords;
msh.MAX =  [max(coords(:, 1)), max(coords(:, 2)), max(coords(:, 3))];
msh.MIN =  [min(coords(:, 1)), min(coords(:, 2)), min(coords(:, 3))];

if lines_idx > 0
    msh.LINES = cell2mat(lines);
end
if lines3_idx > 0
    msh.LINES3 = cell2mat(lines3);
end
if lines4_idx > 0
    msh.LINES4 = cell2mat(lines4);
end
if lines5_idx > 0
    msh.LINES5 = cell2mat(lines5);
end
if lines6_idx > 0
    msh.LINES6 = cell2mat(lines6);
end
if lines7_idx > 0
    msh.LINES7 = cell2mat(lines7);
end
if quads_idx > 0
    msh.QUADS = cell2mat(quads);
end
if quads9_idx > 0
    msh.QUADS9 = cell2mat(quads9);
end
if quads16_idx > 0
    msh.QUADS16 = cell2mat(quads16);
end
if quads25_idx > 0
    msh.QUADS25 = cell2mat(quads25);
end
if quads36_idx > 0
    msh.QUADS36 = cell2mat(quads36);
end
if quads49_idx > 0
    msh.QUADS49 = cell2mat(quads49);
end
if hexas_idx > 0
    msh.HEXAS = cell2mat(hexas);
end
if hexas27_idx > 0
    msh.HEXAS27 = cell2mat(hexas27);
end
if hexas64_idx > 0
    msh.HEXAS64 = cell2mat(hexas64);
end
if hexas125_idx > 0
    msh.HEXAS125 = cell2mat(hexas125);
end
if hexas216_idx > 0
    msh.HEXAS216 = cell2mat(hexas216);
end
if hexas343_idx > 0
    msh.HEXAS343 = cell2mat(hexas343);
end
if points_idx > 0
    msh.PNT = cell2mat(points);
end
% disp('saving mat file');
% if(PhysicalNamesExist)
%     save([output_path,mesh_name,'.mat'],'msh','PhysicalNames','-v7.3');
% else
%     save([output_path,mesh_name,'.mat'],'msh','-v7.3');
% end
%
% end
% ---------------------------------------------------------------------------------------------------

end
