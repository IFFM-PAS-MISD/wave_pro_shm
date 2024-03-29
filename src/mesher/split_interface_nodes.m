function [SemMesh] = split_interface_nodes(SemMesh, Damage)
%
% Split nodes at one or more interfaces; It defines cracks or delaminatios
% Double nodes are introduced and coordinates updated with zero gap
%
% USAGE::
%
%   [SemMesh] = split_interface_nodes(SemMesh,Damage)
%
% Arguments:
%     SemMesh (struct):
%       mesh with difined coordinates, nodes, physical names, etc
%
%     Damage (struct):
%       damage volume and interface elements and interface boundary nodes
%
% Returns:
%     SemMesh (struct):
%       updated mesh with splitted interface nodes
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

% TODO check behaviour for cracks (vertical plane)
for iDamage = 1:Damage.nInterfaces
    [nNodes, ~] = size(SemMesh.nodeCoordinates);
    delamination3Delements = Damage.volume(iDamage).volumeElements;
    delamination2Delements = Damage.interface(iDamage).interfaceElements;
    delamNodes = unique(reshape(SemMesh.element2DNodes(delamination2Delements, :), [], 1));
    outerDelamNodes = unique(reshape(Damage.interfaceBoundary(iDamage).boundaryNodes, [], 1));
    interDelamInterfaceNodes = setdiff(delamNodes, outerDelamNodes);
    % visual check
    % figure
    % plot3(SemMesh.nodeCoordinates(outerDelamNodes,1), SemMesh.nodeCoordinates(outerDelamNodes,2),...
    %      SemMesh.nodeCoordinates(outerDelamNodes,3),'ro');
    % hold on;
    % plot3(SemMesh.nodeCoordinates(interDelamInterfaceNodes,1),...
    %       SemMesh.nodeCoordinates(interDelamInterfaceNodes,2),...
    %       SemMesh.nodeCoordinates(interDelamInterfaceNodes,3),'b.');

    coordsNew = zeros(nNodes + length(interDelamInterfaceNodes), 3);
    coordsNew(1:nNodes, :) = SemMesh.nodeCoordinates;
    SemMesh.nodeCoordinates = coordsNew;
    interDelamInterfaceNodesBelow = interDelamInterfaceNodes;
    interDelamInterfaceNodesAbove = zeros(length(interDelamInterfaceNodes), 1, 'int32');
    for i = 1:length(interDelamInterfaceNodes)
        nodeNum = nNodes + i;
        interDelamInterfaceNodesAbove(i) = int32(nodeNum);
        [elementList, nodeList] = find(SemMesh.element3DNodes(delamination3Delements, :) == interDelamInterfaceNodes(i));
        for j = length(nodeList) / 2 + 1:length(nodeList) % half of the nodes above delamination interface
            % update coordinates
            SemMesh.nodeCoordinates(nodeNum, 1) = ...
               SemMesh.nodeCoordinates(SemMesh.element3DNodes(delamination3Delements(elementList(j)), nodeList(j)), 1);
            SemMesh.nodeCoordinates(nodeNum, 2) = ...
               SemMesh.nodeCoordinates(SemMesh.element3DNodes(delamination3Delements(elementList(j)), nodeList(j)), 2);
            % zero gap
            SemMesh.nodeCoordinates(nodeNum, 3) = ...
                SemMesh.nodeCoordinates(SemMesh.element3DNodes(delamination3Delements(elementList(j)), nodeList(j)), 3);
            % update node number
            SemMesh.element3DNodes(delamination3Delements(elementList(j)), nodeList(j)) = nodeNum;
        end
    end
    SemMesh.Damage(iDamage).interDelamInterfaceNodesBelow = unique(interDelamInterfaceNodesBelow);
    SemMesh.Damage(iDamage).interDelamInterfaceNodesAbove = unique(interDelamInterfaceNodesAbove);
end
% visual check
% iDamage = 1;
% figure
% nodeList = SemMesh.Damage(iDamage).interDelamInterfaceNodesBelow;
% plot3(SemMesh.nodeCoordinates(nodeList,1), SemMesh.nodeCoordinates(nodeList,2),...
%      SemMesh.nodeCoordinates(nodeList,3),'ro');
% hold on;
% nodeList = SemMesh.Damage(iDamage).interDelamInterfaceNodesAbove;
% plot3(SemMesh.nodeCoordinates(nodeList,1),...
%       SemMesh.nodeCoordinates(nodeList,2),...
%       SemMesh.nodeCoordinates(nodeList,3),'b.');

% ---------------------------------------------------------------------------------------------------

end
