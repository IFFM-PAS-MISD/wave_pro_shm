function save_pvd_animation(results_path, nFrames, sampleInterval, timeVector)
%
% Save animation of propagating waves in .pvd Paraview compatible format
% It is an xml format pointing to .vtu files containing all time steps
%
% USAGE::
%
%   save_pvd_animation(results_path,nFrames,sampleInterval,timeVector)
%
% Arguments:
%     results_path (string):
%       folder in which results are stored
%
%     nFrames (integer):
%       number of frames in animation
%
%     sampleInterval (integer):
%       number of samples between consecutive frames
%
%     timeVector (double):
%       time vector corresponding to nSamples
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

file_name = 'wave_pro_shm_results.pvd';
fileID = fopen(fullfile(results_path, 'frames_vtu', filesep, file_name), 'w');
if fileID < 0
    disp('Cannot open file for saving pvd data');
    return
end
% -------------Strings of VTU format---------------------------
StrXmlHeader =      '<?xml version="1.0"?> \n';
StrVTKFileOpen =    '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"> \n';
StrCollectionOpen = '  <Collection> \n';
StrDataSet =        '    <DataSet timestep="%.8e" part="0" file="%s" name="shm_wave_pro_results"/> \n';
StrCollectionClose = '  </Collection> \n';
StrVTKFileClose =   '</VTKFile> \n';
%% write to file
fprintf(fileID, StrXmlHeader);
fprintf(fileID, StrVTKFileOpen);
fprintf(fileID, StrCollectionOpen);
for iFrame = 1:nFrames
    iSample = iFrame * sampleInterval;
    timestep = timeVector(iSample) * 1e6; % [us]
    frameName = ['frame', num2str(iSample, '%07u'), '.vtu'];
    fprintf(fileID, StrDataSet, timestep, frameName);
end
fprintf(fileID, StrCollectionClose);
fprintf(fileID, StrVTKFileClose);
fclose(fileID);
% ---------------------------------------------------------------------------------------------------

end
