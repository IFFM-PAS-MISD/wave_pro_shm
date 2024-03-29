% this is the master script used to run WaveProSHM without GUI

clear all; clc;

addpath(genpath('../src/mesher'));
addpath(genpath('../src/wave_pro_shm/solid3d_engine'));
inputFileName = 'test';

load(['.',filesep,'inputs',filesep,inputFileName]);

% select GPU
gpuDeviceTable
isGPUavailable = true;
isVerbose = true;
gpuDeviceIndex = 3;
selectedGPU = gpuDevice(gpuDeviceIndex);

[isCalculationCanceled,isCalculationFailed] = solid3d_engine(...
                    inputFileName, ...
                    AllData.SemMesh, ...
                    AllData.Materials,...
                    AllData.Excitation, ...
                    AllData.Output,...
                    isGPUavailable, ...
                    isVerbose);

reset(selectedGPU);

rmpath(genpath('../src/mesher'));
rmpath(genpath('../src/wave_pro_shm/solid3d_engine'));