function [Excitation, Output] = retrieve_mesh_transducer_elements(SemMesh, Excitation, Output)
%
% Retrieve from the mesh element numbers and physical names of actuators and sensors
%
%
% USAGE::
%
%   [Excitation,Output] = retrieve_mesh_transducer_elements(SemMesh,Excitation,Output)
%
% Arguments:
%     SemMesh (struct):
%       mesh with difined coordinates, nodes, physical names, etc
%
%     Excitation (struct):
%       tags of volumes/surfaces defining actuators and their electrodes
%
%     Output (struct):
%       tags of volumes/surfaces defining sensors and their electrodes
%
% Returns:
%     Excitation (struct):
%       expanded structure by retrieved elements of actuators
%
%     Output (struct):
%       expanded structure by retrieved elements of sensors
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

% actuators
% extract piezoelectric actuator 3D elements
counter = 0;
for i = 1:length(Excitation.actuatorTagList)
    for j = 1:length(SemMesh.PhysicalVolumes)
        if Excitation.actuatorTagList(i) == SemMesh.PhysicalVolumes(j).physicalTag
            counter = counter + 1;
            Excitation.Actuator(counter).pzt3DElements = SemMesh.PhysicalVolumes(j).elements;
            Excitation.actuatorPhysicalNames{counter} = SemMesh.PhysicalVolumes(j).name;
        end
    end
end
% extract piezoelectric actuator positive electrode 2D elements
counter = 0;
for i = 1:length(Excitation.positiveElectrodeTagList)
    for j = 1:length(SemMesh.PhysicalSurfaces)
        if Excitation.positiveElectrodeTagList(i) == SemMesh.PhysicalSurfaces(j).physicalTag
            counter = counter + 1;
            Excitation.ActuatorPositiveElectrode(counter).pzt2DElements = SemMesh.PhysicalSurfaces(j).elements;
            Excitation.actuatorPositiveElectrodePhysicalNames{counter} = SemMesh.PhysicalSurfaces(j).name;
        end
    end
end
% extract piezoelectric actuator negative electrode 2D elements
counter = 0;
for i = 1:length(Excitation.negativeElectrodeTagList)
    for j = 1:length(SemMesh.PhysicalSurfaces)
        if Excitation.negativeElectrodeTagList(i) == SemMesh.PhysicalSurfaces(j).physicalTag
            counter = counter + 1;
            Excitation.ActuatorNegativeElectrode(counter).pzt2DElements = SemMesh.PhysicalSurfaces(j).elements;
            Excitation.actuatorNegativeElectrodePhysicalNames{counter} = SemMesh.PhysicalSurfaces(j).name;
        end
    end
end
% sensors
% extract piezoelectric actuator 3D elements
counter = 0;
for i = 1:length(Output.sensorTagList)
    for j = 1:length(SemMesh.PhysicalVolumes)
        if Output.sensorTagList(i) == SemMesh.PhysicalVolumes(j).physicalTag
            counter = counter + 1;
            Output.Sensor(counter).pzt3DElements = SemMesh.PhysicalVolumes(j).elements;
            Output.sensorPhysicalNames{counter} = SemMesh.PhysicalVolumes(j).name;
        end
    end
end
% extract piezoelectric sensor positive electrode 2D elements
counter = 0;
for i = 1:length(Output.positiveElectrodeTagList)
    for j = 1:length(SemMesh.PhysicalSurfaces)
        if Output.positiveElectrodeTagList(i) == SemMesh.PhysicalSurfaces(j).physicalTag
            counter = counter + 1;
            Output.SensorPositiveElectrode(counter).pzt2DElements = SemMesh.PhysicalSurfaces(j).elements;
            Output.sensorPositiveElectrodePhysicalNames{counter} = SemMesh.PhysicalSurfaces(j).name;
        end
    end
end
% extract piezoelectric sensor negative electrode 2D elements
counter = 0;
for i = 1:length(Output.negativeElectrodeTagList)
    for j = 1:length(SemMesh.PhysicalSurfaces)
        if Output.negativeElectrodeTagList(i) == SemMesh.PhysicalSurfaces(j).physicalTag
            counter = counter + 1;
            Output.SensorNegativeElectrode(counter).pzt2DElements = SemMesh.PhysicalSurfaces(j).elements;
            Output.sensorNegativeElectrodePhysicalNames{counter} = SemMesh.PhysicalSurfaces(j).name;
        end
    end
end

% ---------------------------------------------------------------------------------------------------

end
