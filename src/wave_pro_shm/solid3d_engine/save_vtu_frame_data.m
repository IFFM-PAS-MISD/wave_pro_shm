function save_vtu_frame_data(model_output_path, iSample, Ux, Uy, Uz, Vx, Vy, Vz, isBinary)
%
% Saves data or rather append data to existing file containing mesh related fields
% in the form of displacements and velocities of a frame of propagating waves in Paraview .vtu 
% compatible format for unstructured grid of arbitrary order Lagrange hexahedron elements
%
% USAGE::
%
%   save_vtu_frame_data(model_output_path,iSample,Ux,Uy,Uz,Vx,Vy,Vz,isBinary)
%
% Arguments:
%     model_output_path (string):
%       folder in which results are stored
%
%     iSample (integer):
%       current sample used for frame number name
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

results_path =  fullfile(model_output_path, filesep, 'frames_vtu', filesep);
if ~exist(results_path, 'dir')
    mkdir(results_path);
end
file_name = ['frame', num2str(iSample, '%07u'), '.vtu'];
fileID = fopen(fullfile(results_path, file_name), 'a', 'ieee-le');
if fileID < 0
    disp('Cannot open file for saving vtu data');
    return
end
% -------------Strings of VTU format---------------------------
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
StrDataArrayClose =      '\n         </DataArray> \n';
StrPointDataClose =       '      </PointData> \n';
StrPieceClose =           '    </Piece> \n';
StrUnstructuredGridClose = '  </UnstructuredGrid> \n';
StrVTKFileClose =         '</VTKFile> \n';
%% write data to file: displacements and velocities
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
