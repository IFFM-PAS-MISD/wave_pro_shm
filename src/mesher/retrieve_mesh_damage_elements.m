function [Damage] = retrieve_mesh_damage_elements(SemMesh, Damage)
%
% Retrieve damage-related element numbers and physical names from the mesh
%
%
% USAGE::
%
%   [Damage] = retrieve_mesh_damage_elements(SemMesh,Damage)
%
% Arguments:
%     SemMesh (struct):
%       mesh with difined coordinates, nodes, physical names, etc
%
%     Damage (struct):
%       tags of volumes and interfaces defining damage
%
% Returns:
%     Damage (struct):
%       expanded structure by retrieved elements
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

% retrieve damage volume 3D elements
Damage.nInterfaces = length(Damage.interfaceTagList);
counter = 0;
for iDamage = 1:Damage.nInterfaces
    for j = 1:length(SemMesh.PhysicalVolumes)
        if Damage.volumeTagList(iDamage) == SemMesh.PhysicalVolumes(j).physicalTag
            counter = counter + 1;
            Damage.volume(counter).volumeElements = SemMesh.PhysicalVolumes(j).elements;
            Damage.damagePhysicalNames{counter} = SemMesh.PhysicalVolumes(j).name;
        end
    end
end

% retrieve damage interface 2D elements
counter = 0;
for iDamage = 1:Damage.nInterfaces
    for j = 1:length(SemMesh.PhysicalSurfaces)
        if Damage.interfaceTagList(iDamage) == SemMesh.PhysicalSurfaces(j).physicalTag
            counter = counter + 1;
            Damage.interface(counter).interfaceElements = SemMesh.PhysicalSurfaces(j).elements;
            Damage.interfacePhysicalNames{counter} = SemMesh.PhysicalSurfaces(j).name;
        end
    end
end
% retrieve damage interface boundary nodes
counter = 0;
for iDamage = 1:Damage.nInterfaces
    for j = 1:length(SemMesh.PhysicalLines)
        if Damage.interfaceBoundaryTagList(iDamage) == SemMesh.PhysicalLines(j).physicalTag
            counter = counter + 1;
            Damage.interfaceBoundary(counter).boundaryNodes = SemMesh.PhysicalLines(j).nodes;
            Damage.boundaryPhysicalNames{counter} = SemMesh.PhysicalLines(j).name;
        end
    end
end

% ---------------------------------------------------------------------------------------------------

end
