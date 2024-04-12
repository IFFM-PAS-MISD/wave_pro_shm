function save_vtu_frame_mesh(model_output_path, iSample, ParaviewMesh, isBinary)
%
% Saves mesh related fields for the frame of propagating waves in Paraview .vtu compatible format
% for unstructured grid of arbitrary order Lagrange hexahedron elements
%
% USAGE::
%
%   save_vtu_frame_mesh(model_output_path,iSample,ParaviewMesh,isBinary)
%
% Arguments:
%     model_output_path (string):
%       folder in which results are stored
%
%     iSample (integer):
%       current sample used for frame number name
%
%     ParaviewMesh (struct):
%       spectral element mesh converted to LAGRANGEHEXAHEDRON
%
%     isBinary (logical):
%       if true - binary file format, if false - ascii format
%
%
% .. Note:: Binary means base64; the implementation is reverse engineered from
%           vtk documentation rather than by using vtkWriter interface
%
% .. Note:: This file uses the base64 encoder from the Apache Commons Codec, 
%           http://commons.apache.org/codec/ and distrubed with MATLAB under
%           the Apache License http://commons.apache.org/license.html
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------
encoder = org.apache.commons.codec.binary.Base64;

LAGRANGEHEXAHEDRON = 72; % element type enumeration: lagrange hexahedron element of arbitrary order
[nElements, nElement3DNodes] = size(ParaviewMesh.element3DNodes);
% number of mesh points (coordinates)
nCoordinatePoints = size(ParaviewMesh.nodeCoordinates, 1);
results_path =  fullfile(model_output_path, filesep, 'frames_vtu', filesep);
if ~exist(results_path, 'dir')
    mkdir(results_path);
end
file_name = ['frame', num2str(iSample, '%07u'), '.vtu'];
fileID = fopen(fullfile(results_path, file_name), 'w', 'ieee-le');
if fileID < 0
    disp('Cannot open file for saving vtu data');
    return
end
% -------------Strings of VTU format---------------------------
StrXmlHeader =           '<?xml version="1.0"?> \n';
StrVTKFileOpen =         '<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32"> \n';
StrUnstructuredGridOpen = '  <UnstructuredGrid> \n';
StrPieceOpen =           '    <Piece NumberOfPoints="%d" NumberOfCells="%d"> \n';
StrPointsOpen =          '      <Points> \n';
if isBinary
    StrDataArrayPoints = '        <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="binary"> \n';
else
    StrDataArrayPoints = '        <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii"> \n';
end
StrDataArrayClose =      '\n         </DataArray> \n';
StrPointsClose =         '      </Points> \n';
StrCellsOpen =           '      <Cells> \n';
if isBinary
    StrDataArrayConnectivity = '        <DataArray type="Int32"  Name="connectivity" format="binary"> \n';
else
    StrDataArrayConnectivity = '        <DataArray type="Int32"  Name="connectivity" format="ascii"> \n';
end
if isBinary
    StrDataArrayOffset =     '         <DataArray type="Int32"  Name="offsets"    format="binary"> \n';
else
    StrDataArrayOffset =     '         <DataArray type="Int32"  Name="offsets"    format="ascii"> \n';
end
if isBinary
    StrDataArrayTypes =      '         <DataArray type="UInt8"  Name="types"      format="binary"> \n';
else
    StrDataArrayTypes =      '         <DataArray type="UInt8"  Name="types"      format="ascii"> \n';
end
StrCellsClose =          '      </Cells> \n';

%% write to file
fprintf(fileID, StrXmlHeader);
fprintf(fileID, StrVTKFileOpen);
fprintf(fileID, StrUnstructuredGridOpen);
fprintf(fileID, StrPieceOpen, nCoordinatePoints, nElements);
fprintf(fileID, StrPointsOpen);
fprintf(fileID, StrDataArrayPoints);
%% coordinates
if isBinary
    % note that transpose of nodeCoordinates is necessary because rows are written first
    dat = reshape(double(ParaviewMesh.nodeCoordinates'), 1, []);
    temp = [encoder.encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')); ...
            encoder.encode(typecast(dat, 'uint8'))];
    fprintf(fileID, '%s', temp);
else
    fprintf(fileID, '%12.8f %12.8f %12.8f \n', ParaviewMesh.nodeCoordinates');
end
fprintf(fileID, StrDataArrayClose);
fprintf(fileID, StrPointsClose);
%% connectivity
fprintf(fileID, StrCellsOpen);
fprintf(fileID, StrDataArrayConnectivity);
% note that transpose of element3DNodes is necessary because rows are written first
if isBinary
    dat = reshape(int32(ParaviewMesh.element3DNodes'), 1, []);
    temp = [encoder.encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')); ...
            encoder.encode(typecast(dat, 'uint8'))];
    fprintf(fileID, '%s', temp);
else
    connectivityFormat = '%d ';
    connectivityFormat = repmat(connectivityFormat, 1, nElement3DNodes);
    connectivityFormat = ['         ', connectivityFormat, '\n'];
    fprintf(fileID, connectivityFormat, ParaviewMesh.element3DNodes');
end
fprintf(fileID, StrDataArrayClose);
%% offset
fprintf(fileID, StrDataArrayOffset);
offset = nElement3DNodes:nElement3DNodes:nElement3DNodes * nElements;
if isBinary
    dat = int32(offset);
    temp = [encoder.encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')); ...
            encoder.encode(typecast(dat, 'uint8'))];
    fprintf(fileID, '%s', temp);
else
    fprintf(fileID, '%d ', offset);
end
fprintf(fileID, StrDataArrayClose);
%% element types
fprintf(fileID, StrDataArrayTypes);
if isBinary
    dat = reshape(uint8(repmat(LAGRANGEHEXAHEDRON, 1, nElements)), 1, []);
    temp = [encoder.encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')); ...
            encoder.encode(typecast(dat, 'uint8'))];
    fprintf(fileID, '%s', temp);
else
    fprintf(fileID, '%d ', repmat(LAGRANGEHEXAHEDRON, 1, nElements));
end
fprintf(fileID, StrDataArrayClose);
fprintf(fileID, StrCellsClose);

fclose(fileID);
% ---------------------------------------------------------------------------------------------------

end
