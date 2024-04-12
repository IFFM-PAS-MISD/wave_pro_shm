function save_vtu_frame(model_output_path, iSample, ParaviewMesh, Ux, Uy, Uz, Vx, Vy, Vz, isBinary)
%
% Saves frame of propagating waves in Paraview .vtu compatible format
% for unstructured grid of arbitrary order Lagrange hexahedron elements
%
% USAGE::
%
%   save_vtu_frame(model_output_path,iSample,ParaviewMesh,Ux,Uy,Uz,Vx,Vy,Vz,isBinary)
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
%     Ui (double):
%       displacement field vector, i=x,y,z
%
%     Vi (double):
%       velocity field vector, i=x,y,z;
%
%     isBinary (logical):
%       if true - binary file format, if false - ascii format
%
%
% .. Note:: Binary means base64; the implementation is reverse engineered from
%           vtk documentation rather than by using vtkWriter interface
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

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
StrPointDataOpen =       '      <PointData> \n';
if isBinary
    StrDataArrayDisplX =  '         <DataArray type="Float64" Name="DisplacementsX" NumberOfComponents="1" format="binary"> \n';
    StrDataArrayDisplY =  '         <DataArray type="Float64" Name="DisplacementsY" NumberOfComponents="1" format="binary"> \n';
    StrDataArrayDisplZ =  '         <DataArray type="Float64" Name="DisplacementsZ" NumberOfComponents="1" format="binary"> \n';
    StrDataArrayVelocX =  '         <DataArray type="Float64" Name="VelocitiesX" NumberOfComponents="1" format="binary"> \n';
    StrDataArrayVelocY =  '         <DataArray type="Float64" Name="VelocitiesY" NumberOfComponents="1" format="binary"> \n';
    StrDataArrayVelocZ =  '         <DataArray type="Float64" Name="VelocitiesZ" NumberOfComponents="1" format="binary"> \n';
else
    StrDataArrayDisplX =  '         <DataArray type="Float64" Name="DisplacementsX" NumberOfComponents="1" format="ascii"> \n';
    StrDataArrayDisplY =  '         <DataArray type="Float64" Name="DisplacementsY" NumberOfComponents="1" format="ascii"> \n';
    StrDataArrayDisplZ =  '         <DataArray type="Float64" Name="DisplacementsZ" NumberOfComponents="1" format="ascii"> \n';
    StrDataArrayVelocX =  '         <DataArray type="Float64" Name="VelocitiesX" NumberOfComponents="1" format="ascii"> \n';
    StrDataArrayVelocY =  '         <DataArray type="Float64" Name="VelocitiesY" NumberOfComponents="1" format="ascii"> \n';
    StrDataArrayVelocZ =  '         <DataArray type="Float64" Name="VelocitiesZ" NumberOfComponents="1" format="ascii"> \n';
end
StrPointDataClose =       '      </PointData> \n';
StrPieceClose =           '    </Piece> \n';
StrUnstructuredGridClose = '  </UnstructuredGrid> \n';
StrVTKFileClose =         '</VTKFile> \n';
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
    temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
            matlab.net.base64encode(typecast(dat, 'uint8'))];
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
    temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
            matlab.net.base64encode(typecast(dat, 'uint8'))];
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
    temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
            matlab.net.base64encode(typecast(dat, 'uint8'))];
    fprintf(fileID, '%s', temp);
else
    fprintf(fileID, '%d ', offset);
end
fprintf(fileID, StrDataArrayClose);
%% element types
fprintf(fileID, StrDataArrayTypes);
if isBinary
    dat = reshape(uint8(repmat(LAGRANGEHEXAHEDRON, 1, nElements)), 1, []);
    temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
            matlab.net.base64encode(typecast(dat, 'uint8'))];
    fprintf(fileID, '%s', temp);
else
    fprintf(fileID, '%d ', repmat(LAGRANGEHEXAHEDRON, 1, nElements));
end
fprintf(fileID, StrDataArrayClose);
fprintf(fileID, StrCellsClose);
%% data: displacements and velocities
fprintf(fileID, StrPointDataOpen);
% displacements X
if ~isempty(Ux)
    fprintf(fileID, StrDataArrayDisplX);
    if isBinary
        Ux(isnan(Ux)) = 0;
        dat = double(Ux);
        temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
                matlab.net.base64encode(typecast(dat, 'uint8'))];
        fprintf(fileID, '%s', temp);
    else
        fprintf(fileID, '         %.15e \n', Ux);
    end
    fprintf(fileID, StrDataArrayClose);
end
% displacements Y
if ~isempty(Uy)
    fprintf(fileID, StrDataArrayDisplY);
    if isBinary
        Uy(isnan(Uy)) = 0;
        dat = double(Uy);
        temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
                matlab.net.base64encode(typecast(dat, 'uint8'))];
        fprintf(fileID, '%s', temp);
    else
        fprintf(fileID, '         %.15e \n', Uy);
    end
    fprintf(fileID, StrDataArrayClose);
end
% displacements Z
if ~isempty(Uz)
    fprintf(fileID, StrDataArrayDisplZ);
    if isBinary
        Uz(isnan(Uz)) = 0;
        dat = double(Uz);
        temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
                matlab.net.base64encode(typecast(dat, 'uint8'))];
        fprintf(fileID, '%s', temp);
    else
        fprintf(fileID, '         %.15e \n', Uz);
    end
    fprintf(fileID, StrDataArrayClose);
end
% velocities X
if ~isempty(Vx)
    fprintf(fileID, StrDataArrayVelocX);
    if isBinary
        Vx(isnan(Vx)) = 0;
        dat = double(Vx);
        temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
                matlab.net.base64encode(typecast(dat, 'uint8'))];
        fprintf(fileID, '%s', temp);
    else
        fprintf(fileID, '         %.15e \n', Vx);
    end
    fprintf(fileID, StrDataArrayClose);
end
% velocities Y
if ~isempty(Vy)
    fprintf(fileID, StrDataArrayVelocY);
    if isBinary
        Vy(isnan(Vy)) = 0;
        dat = double(Vy);
        temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
                matlab.net.base64encode(typecast(dat, 'uint8'))];
        fprintf(fileID, '%s', temp);
        fprintf(fileID, '%s', temp);
    else
        fprintf(fileID, '         %.15e \n', Vy);
    end
    fprintf(fileID, StrDataArrayClose);
end
% velocities Z
if ~isempty(Vz)
    fprintf(fileID, StrDataArrayVelocZ);
    if isBinary
        Vz(isnan(Vz)) = 0;
        dat = double(Vz);
        temp = [matlab.net.base64encode(typecast(uint32(numel(dat) * numel(typecast(dat(1), 'uint8'))), 'uint8')), ...
                matlab.net.base64encode(typecast(dat, 'uint8'))];
        fprintf(fileID, '%s', temp);
    else
        fprintf(fileID, '         %.15e \n', Vz);
    end
    fprintf(fileID, StrDataArrayClose);
end
fprintf(fileID, StrPointDataClose);
fprintf(fileID, StrPieceClose);
fprintf(fileID, StrUnstructuredGridClose);
fprintf(fileID, StrVTKFileClose);
fclose(fileID);
% ---------------------------------------------------------------------------------------------------

end
