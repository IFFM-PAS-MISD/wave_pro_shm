function [isCalculationCanceled, isCalculationFailed] = solid3d_engine(...
          inputFileName, SemMesh, Materials, Excitation, Output, ...
          isGPUavailable, isVerbose, SignalPreview, calculationProgress)
%
% Engine of solid 3D spectral element method
% It computes elastic wave propagation in isotropic or orthotropic structural elements
% Current implementation is limited to materials aligned to x-y plane
%
% USAGE::
%
%   [isCalculationCanceled,isCalculationFailed] = solid3d_engine(...
%    inputFileName, SemMesh, Materials,Excitation, Output,...
%    isGPUavailable, isVerbose, SignalPreview, calculationProgress)
%
% Arguments:
%     inputFileName (string):
%       the name of input file; the output folder inherits this name
%
%     SemMesh (struct):
%       spectral element mesh
%
%     Materials (struct):
%       material properties and their assignment to volumes
%
%     Excitation (struct):
%       excitation signals and their assignment to actuators
%
%     Output (struct):
%       sensors, response points and output related infos
%
%     isGPUavailable (logical):
%       flag for gpu availability, if true - computations are on gpu
%                                  if false - computations are on cpu
%
%     isVerbose (logical):
%       flag for spitting additional information during calculations
%
%     SignalPreview (integer):
%       handle to the app signal preview axes for updating signal plot
%
%     calculationProgress (integer):
%       handle to the app progress dialog box
%
% .. Note::
%    Function can take 7 arguments instead of 9 for running without GUI
%
% Returns:
%     isCalculationCanceled (logical):
%       flag for checking if user cancelled calculations in gui
%
%     isCalculationFailed (logical):
%       flag is true if something goes wrong
%       e.g. integration error or problem with saving file
%
%
% .. seealso:: Functions :func:`gll`,
%                        :func:`vandermonde`,
%                        :func:`shape3d_prim_coeff`,
%                        :func:`sparse_block_diagonal_shape_derivatives`,
%                        :func:`inv_jacp`,
%                        :func:`det_jacp`,
%                        :func:`local_global_nodes_map2dofs`,
%                        :func:`excitation_signals_at_actuators`,
%                        :func:`retrieve_mesh_transducer_elements`,
%                        :func:`pzt_coupling_matrices`,
%                        :func:`dofs3d`,
%                        :func:`dofs3dfifi`,
%                        :func:`sem2paraview_mesh`,
%                        :func:`strains_spec_p`,
%                        :func:`stresses_spec_p`,
%                        :func:`save_vtu_frame`,
%                        :func:`save_pvd_animation`
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

if(nargin == 7)
    SignalPreview = [];
    calculationProgress = [];
end

isCalculationCanceled = false;
isCalculationFailed = false;
model_output_path = fullfile('.', filesep, 'outputs', filesep, inputFileName, filesep);
if ~exist(model_output_path, 'dir')
    mkdir(model_output_path);
end

[nElements, nElement3DNodes] = size(SemMesh.element3DNodes);
x1 = SemMesh.nodeCoordinates(SemMesh.element3DNodes(1, 1), 1);
x2 = SemMesh.nodeCoordinates(SemMesh.element3DNodes(1, 2), 1);
y1 = SemMesh.nodeCoordinates(SemMesh.element3DNodes(1, 1), 2);
y2 = SemMesh.nodeCoordinates(SemMesh.element3DNodes(1, 2), 2);
z1 = SemMesh.nodeCoordinates(SemMesh.element3DNodes(1, 1), 3);
z2 = SemMesh.nodeCoordinates(SemMesh.element3DNodes(1, 2), 3);
characteristicDistance = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
%% assign materials - elasticity matrix
% _P at the end means that the matrices and vectors are stored by disconected elements
% for parallel computation
if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Material assignment';
    else
        disp('Material assignment');
    end
end

if isGPUavailable  % allocate matrices directly on gpu
    D11_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D12_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D13_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D14_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D22_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D23_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D24_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D33_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D34_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D44_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D55_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D56_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    D66_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');

    rho_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    dampingCoeffX_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    dampingCoeffY_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    dampingCoeffZ_P = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
else
    D11_P = zeros(nElements * nElement3DNodes, 1);
    D12_P = zeros(nElements * nElement3DNodes, 1);
    D13_P = zeros(nElements * nElement3DNodes, 1);
    D14_P = zeros(nElements * nElement3DNodes, 1);
    D22_P = zeros(nElements * nElement3DNodes, 1);
    D23_P = zeros(nElements * nElement3DNodes, 1);
    D24_P = zeros(nElements * nElement3DNodes, 1);
    D33_P = zeros(nElements * nElement3DNodes, 1);
    D34_P = zeros(nElements * nElement3DNodes, 1);
    D44_P = zeros(nElements * nElement3DNodes, 1);
    D55_P = zeros(nElements * nElement3DNodes, 1);
    D56_P = zeros(nElements * nElement3DNodes, 1);
    D66_P = zeros(nElements * nElement3DNodes, 1);

    rho_P = zeros(nElements * nElement3DNodes, 1);
    dampingCoeffX_P = zeros(nElements * nElement3DNodes, 1);
    dampingCoeffY_P = zeros(nElements * nElement3DNodes, 1);
    dampingCoeffZ_P = zeros(nElements * nElement3DNodes, 1);
end

for iVolume = 1:length(SemMesh.PhysicalVolumes)
    % Check for Cancel button press
    if(~isempty(calculationProgress))
        if calculationProgress.CancelRequested
            isCalculationCanceled = true;
            return
        end
    end
    iMaterial = Materials.assignment(iVolume).index;

    % isotropic material or orthotropic material
    if Materials.assignment(iVolume).category == 1 || ...
       Materials.assignment(iVolume).category == 2
        elementList = SemMesh.PhysicalVolumes(iVolume).elements;
        for j = 1:length(elementList)
            elementNo = elementList(j);
            n1 = (elementNo - 1) * nElement3DNodes + 1;
            n2 = n1 + nElement3DNodes - 1;
            D11_P(n1:n2, 1) = Materials.elasticityMatrix(1, 1, iMaterial);
            D12_P(n1:n2, 1) = Materials.elasticityMatrix(1, 2, iMaterial);
            D13_P(n1:n2, 1) = Materials.elasticityMatrix(1, 3, iMaterial);
            D14_P(n1:n2, 1) = Materials.elasticityMatrix(1, 4, iMaterial);
            D22_P(n1:n2, 1) = Materials.elasticityMatrix(2, 2, iMaterial);
            D23_P(n1:n2, 1) = Materials.elasticityMatrix(2, 3, iMaterial);
            D24_P(n1:n2, 1) = Materials.elasticityMatrix(2, 4, iMaterial);
            D33_P(n1:n2, 1) = Materials.elasticityMatrix(3, 3, iMaterial);
            D34_P(n1:n2, 1) = Materials.elasticityMatrix(3, 4, iMaterial);
            D44_P(n1:n2, 1) = Materials.elasticityMatrix(4, 4, iMaterial);
            D55_P(n1:n2, 1) = Materials.elasticityMatrix(5, 5, iMaterial);
            D56_P(n1:n2, 1) = Materials.elasticityMatrix(5, 6, iMaterial);
            D66_P(n1:n2, 1) = Materials.elasticityMatrix(6, 6, iMaterial);
            rho_P(n1:n2, 1) = Materials.density(1, iMaterial);
            dampingCoeffX_P(n1:n2, 1) = Materials.damping(1, iMaterial);
            dampingCoeffY_P(n1:n2, 1) = Materials.damping(2, iMaterial);
            dampingCoeffZ_P(n1:n2, 1) = Materials.damping(3, iMaterial);
        end
    elseif Materials.assignment(iVolume).category == 3  % piezoelectric material
        elementList = SemMesh.PhysicalVolumes(iVolume).elements;
        for j = 1:length(elementList)
            elementNo = elementList(j);
            n1 = (elementNo - 1) * nElement3DNodes + 1;
            n2 = n1 + nElement3DNodes - 1;
            D11_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(1, 1, iMaterial);
            D12_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(1, 2, iMaterial);
            D13_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(1, 3, iMaterial);
            D14_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(1, 4, iMaterial);
            D22_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(2, 2, iMaterial);
            D23_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(2, 3, iMaterial);
            D24_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(2, 4, iMaterial);
            D33_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(3, 3, iMaterial);
            D34_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(3, 4, iMaterial);
            D44_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(4, 4, iMaterial);
            D55_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(5, 5, iMaterial);
            D56_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(5, 6, iMaterial);
            D66_P(n1:n2, 1) = Materials.elasticityMatrixPiezo(6, 6, iMaterial);
            rho_P(n1:n2, 1) = Materials.densityPiezo(1, iMaterial);
            dampingCoeffX_P(n1:n2, 1) = Materials.dampingPiezo(1, iMaterial);
            dampingCoeffY_P(n1:n2, 1) = Materials.dampingPiezo(2, iMaterial);
            dampingCoeffZ_P(n1:n2, 1) = Materials.dampingPiezo(3, iMaterial);
        end

    end
end
if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Material assignment. Done. ';
    else
        disp('Material assignment. Done. ');
    end
end
% number of nodes along local x axis of element (ksi-axis)
nElementNodesKsi = SemMesh.elementOrder + 1;
% number of nodes along local y axis of element (eta-axis)
nElementNodesEta = SemMesh.elementOrder + 1;
% number of nodes along local z axis of element (dzeta-axis)
nElementNodesDzeta = SemMesh.elementOrder + 1;
[ksi, wKsi] = gll(nElementNodesKsi); % weights and nodes distribution
[eta, wEta] = gll(nElementNodesEta);
[dzeta, wDzeta] = gll(nElementNodesDzeta);
% vandermonde approach
vandermondeMatrixKsi = vandermonde(ksi);
vandermondeMatrixEta = vandermonde(eta);
vandermondeMatrixDzeta = vandermonde(dzeta);

inverseVandermondeMatrixKsi = inv(vandermondeMatrixKsi);
inverseVandermondeMatrixEta = inv(vandermondeMatrixEta);
inverseVandermondeMatrixDzeta = inv(vandermondeMatrixDzeta);

% if(isVerbose)
%     disp('local derivatives (coefficients for each shape function at each node)');
% end

if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Local derivatives';
    else
        disp('Local derivatives');
    end
end

% create sparse block-diagonal matrix of shape function derivatives
% _e at the end means that this is for a single element
[shapePrimKsi_e, shapePrimEta_e, shapePrimDzeta_e] = shape3d_prim_coeff(nElementNodesKsi, ...
                                                                        nElementNodesEta, nElementNodesDzeta, inverseVandermondeMatrixKsi, ...
                                                                        inverseVandermondeMatrixEta, inverseVandermondeMatrixDzeta, ksi', eta', dzeta');
% abbreviate naming: Nksi is shape function derivative in respect to ksi
[Nksi_P, Neta_P, Ndzeta_P] = sparse_block_diagonal_shape_derivatives(shapePrimKsi_e, shapePrimEta_e, ...
                                                                     shapePrimDzeta_e, nElements, nElement3DNodes);
% global coordinates
X = zeros(nElements * nElement3DNodes, 1);
Y = zeros(nElements * nElement3DNodes, 1);
Z = zeros(nElements * nElement3DNodes, 1);

for iElement = 1:nElements
    n1 = (iElement - 1) * nElement3DNodes + 1;
    n2 = n1 + nElement3DNodes - 1;
    X(n1:n2) = SemMesh.nodeCoordinates(SemMesh.element3DNodes(iElement, :), 1);
    Y(n1:n2) = SemMesh.nodeCoordinates(SemMesh.element3DNodes(iElement, :), 2);
    Z(n1:n2) = SemMesh.nodeCoordinates(SemMesh.element3DNodes(iElement, :), 3);
end
if isGPUavailable  % push to gpu
    try
        Nksi_P = gpuArray(Nksi_P);
    catch
        if isVerbose
            disp('Cannot push data to GPU, please check if GPU is available');
        end
        isCalculationFailed = true;
        return
    end
    Neta_P = gpuArray(Neta_P);
    Ndzeta_P = gpuArray(Ndzeta_P);
    X = gpuArray(X);
    Y = gpuArray(Y);
    Z = gpuArray(Z);
end

% NksiT means transpose
NksiT_P = Nksi_P';
NetaT_P = Neta_P';
NdzetaT_P = Ndzeta_P';

if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Jacobians';
    else
        disp('Jacobians');
    end
end
J11_P = Nksi_P * X;
J21_P = Nksi_P * Y;
J31_P = Nksi_P * Z;
J12_P = Neta_P * X;
J22_P = Neta_P * Y;
J32_P = Neta_P * Z;
J13_P = Ndzeta_P * X;
J23_P = Ndzeta_P * Y;
J33_P = Ndzeta_P * Z;

% if(isVerbose)
%     disp('determinant and inverse Jacobians')
% end
if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Determinants and inverse Jacobians';
    else
        disp('Determinants and inverse Jacobians');
    end
end
% /* determinant and inverse Jacobians  */
[invJ11_P, invJ12_P, invJ13_P, invJ21_P, invJ22_P, invJ23_P, invJ31_P, invJ32_P, invJ33_P] = ...
    inv_jacp(J11_P, J12_P, J13_P, J21_P, J22_P, J23_P, J31_P, J32_P, J33_P);

detJ_P = det_jacp(J11_P, J12_P, J13_P, J21_P, J22_P, J23_P, J31_P, J32_P, J33_P);

negativeJnodeList = find(detJ_P < 0);
negativeJelementList = unique(floor(negativeJnodeList / nElement3DNodes) + 1);

if ~isempty(negativeJelementList)  % check for negative Jacobian
    if isVerbose
        if(~isempty(calculationProgress))
            calculationProgress.Value = 0;
            calculationProgress.Message = 'Error: negative Jacobians detected';
        else
            disp('Error: negative Jacobians detected');
        end
        return
    end
end

% weights at quadrature points
counter = 0;
www_e = zeros(nElement3DNodes, 1);
for kz = 1:nElementNodesDzeta % dzeta
    for ky = 1:nElementNodesEta % eta
        for kx = 1:nElementNodesKsi % ksi
            counter = counter + 1;
            www_e(counter, 1) = wKsi(kx) * wEta(ky) * wDzeta(kz);
        end
    end
end
if isGPUavailable  % push to gpu
    www_e = gpuArray(www_e);
end
%
WWW_P = repmat(www_e, nElements, 1);
WWWDetJ_P = WWW_P .* detJ_P;

clear WWW www detJ_P;
clear J11_P J12_P J13_P J21_P J22_P J23_P J31_P J32_P J33_P;

% mass matrix nodal values (vector because of diagonality)
if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Mass matrix';
    else
        disp('Mass matrix');
    end
end
massVector_P = rho_P .* WWWDetJ_P;

if isGPUavailable  % allocate matrices directly on gpu
    massMatrixGlobal = zeros(SemMesh.nDofs, 1, 'double', 'gpuArray');
    dampingMatrixGlobal = zeros(SemMesh.nDofs, 1, 'double', 'gpuArray');
    massVectorByDof_P = zeros(3 * nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    dampingVectorByDof_P = zeros(3 * nElements * nElement3DNodes, 1, 'double', 'gpuArray');
else    % allocate matrices on cpu
    massMatrixGlobal = zeros(SemMesh.nDofs, 1);
    dampingMatrixGlobal = zeros(SemMesh.nDofs, 1);
    massVectorByDof_P = zeros(3 * nElements * nElement3DNodes, 1);
    dampingVectorByDof_P = zeros(3 * nElements * nElement3DNodes, 1);
end

massVectorByDof_P(1:3:end, :) = massVector_P;
massVectorByDof_P(2:3:end, :) = massVector_P;
massVectorByDof_P(3:3:end, :) = massVector_P;

dampingVectorByDof_P(1:3:end, :) = dampingCoeffX_P .* massVector_P;
dampingVectorByDof_P(2:3:end, :) = dampingCoeffY_P .* massVector_P;
dampingVectorByDof_P(3:3:end, :) = dampingCoeffZ_P .* massVector_P;

% assembly mass and damping matrices
[I_G_dofs, I_L_dofs] = local_global_nodes_map2dofs(SemMesh.I_G, SemMesh.I_L, isGPUavailable);

for ibasket = 1:SemMesh.nBaskets
    massMatrixGlobal(I_G_dofs(:, ibasket)) = massMatrixGlobal(I_G_dofs(:, ibasket)) + massVectorByDof_P(I_L_dofs(:, ibasket));
    dampingMatrixGlobal(I_G_dofs(:, ibasket)) = dampingMatrixGlobal(I_G_dofs(:, ibasket)) + dampingVectorByDof_P(I_L_dofs(:, ibasket));
end

clear massVector_P massVectorByDof_P rho_P dampingVector_P dampingVectorByDof_P;
clear dampingCoeffX_P dampingCoeffY_P dampingCoeffZ_P;

% degrees of freedom to nodes translation for storage by disconnected elemetns
indexMapX = zeros(nElements * nElement3DNodes, 1);
indexMapY = zeros(nElements * nElement3DNodes, 1);
indexMapZ = zeros(nElements * nElement3DNodes, 1);
for iElement = 1:nElements
    n1 = (iElement - 1) * nElement3DNodes + 1;
    n2 = n1 + nElement3DNodes - 1;
    indexMapX(n1:n2, 1) = 3 * SemMesh.element3DNodes(iElement, :) - 2;
    indexMapY(n1:n2, 1) = 3 * SemMesh.element3DNodes(iElement, :) - 1;
    indexMapZ(n1:n2, 1) = 3 * SemMesh.element3DNodes(iElement, :);
end
if isGPUavailable  % push to gpu
    indexMapX = gpuArray(uint32(indexMapX));
    indexMapY = gpuArray(uint32(indexMapY));
    indexMapZ = gpuArray(uint32(indexMapZ));
end

%% piezoelectric transducers / signals
if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'PZT actuators';
    else
        disp('PZT actuators');
    end
end
% Excitation signals for all actuators
Excitation = excitation_signals_at_actuators(Excitation);
% global vector of forces induced by the actuator
% forcesActuators = zeros(SemMesh.nDofs,Excitation.nActuators);

% retrieve elements of actuators and sensors from mesh
[Excitation, Output] = retrieve_mesh_transducer_elements(SemMesh, Excitation, Output);

% TODO improve the code - it is slow here!
% TODO add electrode surfaces to mesh PhysicalNames or implement better method for extraction
% of grounded and top electrodes - current version will not work when multiple elements are stack
% on top of each other; it is only for one layer of elements per actuator/sensor

% global matrices
% notations
% fi - electric potential
% u - displacements
% for electric boundary conditions:
% 0 - grounded nodes or dofs
% i - internal nodes
% n - nodes or dofs related to positive electrode where voltage is applied

Kfifi = sparse(SemMesh.nNodes, SemMesh.nNodes);
Kfiu  = sparse(SemMesh.nNodes, SemMesh.nDofs);

for iActuator = 1:Excitation.nActuators
    % Check for Cancel button press
    if(~isempty(calculationProgress))
        if calculationProgress.CancelRequested
            isCalculationCanceled = true;
            return
        end
    end
    iVolume = Excitation.actuatorTagList(iActuator);
    iPiezoMat = Materials.assignment(iVolume).index;
    if isVerbose
        if(~isempty(calculationProgress))
            calculationProgress.Value = 0;
            calculationProgress.Message = ['Actuator ', Excitation.actuatorPhysicalNames{iActuator}];
        else
            disp(['Actuator ', Excitation.actuatorPhysicalNames{iActuator}]);
        end
    end
    elementList = Excitation.Actuator(iActuator).pzt3DElements;
    nodespzt = SemMesh.element3DNodes(elementList, :);
    % electric boundary conditions
    % bottom surface
    bottomSurfElements = Excitation.ActuatorNegativeElectrode(iActuator).pzt2DElements;
    EBC0 =  reshape(SemMesh.element2DNodes(bottomSurfElements, :), [], 1);
    EBC0 = unique(EBC0);
    % top surface
    topSurfElements = Excitation.ActuatorPositiveElectrode(iActuator).pzt2DElements;
    EBCN = reshape(SemMesh.element2DNodes(topSurfElements, :), [], 1);
    EBCN = unique(EBCN);
    % internal nodes
    EBCI = setdiff(reshape(nodespzt, [], 1), union(EBC0, EBCN));
    EBCI = unique(EBCI);

    BC03 = zeros(3 * length(EBC0), 1);
    BC03(1:3:end, 1) = 3 * EBC0 - 2;
    BC03(2:3:end, 1) = 3 * EBC0 - 1;
    BC03(3:3:end, 1) = 3 * EBC0 - 0;

    BCI3 = zeros(3 * length(EBCI), 1);
    BCI3(1:3:end, 1) = 3 * EBCI - 2;
    BCI3(2:3:end, 1) = 3 * EBCI - 1;
    BCI3(3:3:end, 1) = 3 * EBCI - 0;

    BCN3 = zeros(3 * length(EBCN), 1);
    BCN3(1:3:end, 1) = 3 * EBCN - 2;
    BCN3(2:3:end, 1) = 3 * EBCN - 1;
    BCN3(3:3:end, 1) = 3 * EBCN - 0;

    element3DDofspzt = zeros(length(elementList), 3 * nElement3DNodes);

    element3DDofspzt(:, 1:3:end) = 3 * nodespzt(:, :) - 2;
    element3DDofspzt(:, 2:3:end) = 3 * nodespzt(:, :) - 1;
    element3DDofspzt(:, 3:3:end) = 3 * nodespzt(:, :) - 0;

    cpzt = 0;
    for elementNo = elementList   % local matrices for each element in pzt
        cpzt = cpzt + 1;
        nodeList = SemMesh.element3DNodes(elementNo, :);
        [kufi, kfifi] = pzt_coupling_matrices(Materials.pztCoupling(:, :, iPiezoMat), ...
                                              Materials.pztPermittivity(:, :, iPiezoMat), ...
                                              SemMesh.nodeCoordinates(nodeList, 1), SemMesh.nodeCoordinates(nodeList, 2), ...
                                              SemMesh.nodeCoordinates(nodeList, 3), inverseVandermondeMatrixKsi, ...
                                              inverseVandermondeMatrixEta, inverseVandermondeMatrixDzeta, ...
                                              ksi, eta, dzeta, wKsi, wEta, wDzeta);

        kfiu = sparse(kufi');
        kfifi = sparse(kfifi);

        [iout] = dofs3d(element3DDofspzt, cpzt, SemMesh.nDofs, nElement3DNodes);
        [iout2] = dofs3dfifi(nodespzt, cpzt, SemMesh.nNodes, nElement3DNodes);

        Kfiu = Kfiu + iout2' * kfiu * iout; % assembly
        Kfifi = Kfifi + iout2' * kfifi * iout2; % assembly
    end

    % KfiuA = [Kfiu(EBCI,BC03),Kfiu(EBCI,BCI3),Kfiu(EBCI,BCN3)];% actuator BC
    % KfifiA = Kfifi(EBCI,EBCI); % actuator
    % Kfifiin = Kfifi(EBCI,EBCN);
    % V1 = sparse(SemMesh.nNodes,1);
    % if(Excitation.assignment(iActuator).category == 1) % Hann signal
    %     iHannSignal = Excitation.assignment(iActuator).index;
    %     V1(EBCN,1) = Excitation.Hann(iHannSignal).peakVoltage; % top
    % end
    % V1(EBC0,1) = 0; % bottom
    % fiAi = -KfifiA \ (Kfifiin * V1(EBCN)); % induced potential
    % FiA = sparse(SemMesh.nNodes,1);
    %
    % if(Excitation.assignment(iActuator).category == 1) % Hann signal
    %     iHannSignal = Excitation.assignment(iActuator).index;
    %     FiA(EBCN) = Excitation.Hann(iHannSignal).peakVoltage;
    % end
    % FiA(EBCI) = fiAi;
    % FaL = Kfiu' * FiA; % local forces in pzt actuator of index iActuator

    % BN = [BC03;BCI3;BCN3];
    % forcesActuators(BN,iActuator) = FaL(BN); % global forces from pzt actuator of index iActuator
    PztActuator(iActuator).EBCI = EBCI;
    PztActuator(iActuator).EBCIN = [EBCI; EBCN];
    PztActuator(iActuator).EBCN = EBCN;
    PztActuator(iActuator).EBC0 = EBC0;
    PztActuator(iActuator).BC03 = BC03;
    PztActuator(iActuator).BCI3 = BCI3;
    PztActuator(iActuator).BCN3 = BCN3;
end
actuatorsEBC0 = []; % grounded electrode nodes
actuatorsEBCI = []; % internal nodes
actuatorsEBCN = []; % top electrode nodes
actuatorsEBCIN = []; % internal nodes + top electrode nodes
actuatorsBC03 = []; % grounded electrode degrees of freedom
actuatorsBCI3 = []; % internal element degrees of freedom
actuatorsBCN3 = []; % top electrode degrees of freedom

for iActuator = 1:Excitation.nActuators
    actuatorsEBC0 = [actuatorsEBC0; PztActuator(iActuator).EBC0];
    actuatorsEBCI = [actuatorsEBCI; PztActuator(iActuator).EBCI];
    actuatorsEBCN = [actuatorsEBCN; PztActuator(iActuator).EBCN];
    actuatorsBC03 = [actuatorsBC03; PztActuator(iActuator).BC03];
    actuatorsBCI3 = [actuatorsBCI3; PztActuator(iActuator).BCI3];
    actuatorsBCN3 = [actuatorsBCN3; PztActuator(iActuator).BCN3];

end
actuatorsBC = [actuatorsBC03; actuatorsBCI3; actuatorsBCN3];
actuatorsEBC = [actuatorsEBC0; actuatorsEBCI; actuatorsEBCN];
KfifiAi = Kfifi(actuatorsEBCI, actuatorsEBCI);
invKfifiAi = inv(KfifiAi);
invKfifiAi = full(invKfifiAi);
KfiuAi = [Kfiu(actuatorsEBCI, actuatorsBC03), Kfiu(actuatorsEBCI, actuatorsBCI3), Kfiu(actuatorsEBCI, actuatorsBCN3)]; % all actuators BC
KfiuA = [Kfiu(actuatorsEBC0, actuatorsBC03), Kfiu(actuatorsEBC0, actuatorsBCI3), Kfiu(actuatorsEBC0, actuatorsBCN3)
         Kfiu(actuatorsEBCI, actuatorsBC03), Kfiu(actuatorsEBCI, actuatorsBCI3), Kfiu(actuatorsEBCI, actuatorsBCN3)
         Kfiu(actuatorsEBCN, actuatorsBC03), Kfiu(actuatorsEBCN, actuatorsBCI3), Kfiu(actuatorsEBCN, actuatorsBCN3)];
KfifiAin = Kfifi(actuatorsEBCI, actuatorsEBCN);
% clear FaL
if isGPUavailable  % push to gpu
    % forcesActuators = gpuArray(forcesActuators);
    actuatorsBC = gpuArray(actuatorsBC);
    KfifiAi = gpuArray(KfifiAi);
    invKfifiAi = gpuArray(invKfifiAi);
    KfiuAi = gpuArray(KfiuAi);
    KfifiAin = gpuArray(KfifiAin);
end
%% prepare matrices for sensing
if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'PZT sensors';
    else
        disp('PZT sensors');
    end
end

% global matrices
Kfifi = sparse(SemMesh.nNodes, SemMesh.nNodes);
Kfiu = sparse(SemMesh.nNodes, SemMesh.nDofs);
for iSensor = 1:Output.nSensors
    % Check for Cancel button press
    if(~isempty(calculationProgress))
        if calculationProgress.CancelRequested
            isCalculationCanceled = true;
            return
        end
    end
    iVolume = Output.sensorTagList(iSensor);
    iPiezoMat = Materials.assignment(iVolume).index;
    if isVerbose
        if(~isempty(calculationProgress))
            calculationProgress.Value = 0;
            calculationProgress.Message = ['Sensor ', Output.sensorPhysicalNames{iSensor}];
        else
            disp(['Sensor ', Output.sensorPhysicalNames{iSensor}]);
        end
    end
    elementList = Output.Sensor(iSensor).pzt3DElements;
    nodespzt = SemMesh.element3DNodes(elementList, :);
    % electric boundary conditions
    % bottom surface
    bottomSurfElements = Output.SensorNegativeElectrode(iSensor).pzt2DElements;
    EBC0 =  reshape(SemMesh.element2DNodes(bottomSurfElements, :), [], 1);
    EBC0 = unique(EBC0);
    % top surface
    topSurfElements = Output.SensorPositiveElectrode(iSensor).pzt2DElements;
    EBCN = reshape(SemMesh.element2DNodes(topSurfElements, :), [], 1);
    EBCN = unique(EBCN);
    % internal nodes
    EBCI = setdiff(reshape(nodespzt, [], 1), union(EBC0, EBCN));
    EBCI = unique(EBCI);

    BC03 = zeros(3 * length(EBC0), 1);
    BC03(1:3:end, 1) = 3 * EBC0 - 2;
    BC03(2:3:end, 1) = 3 * EBC0 - 1;
    BC03(3:3:end, 1) = 3 * EBC0 - 0;

    BCI3 = zeros(3 * length(EBCI), 1);
    BCI3(1:3:end, 1) = 3 * EBCI - 2;
    BCI3(2:3:end, 1) = 3 * EBCI - 1;
    BCI3(3:3:end, 1) = 3 * EBCI - 0;

    BCN3 = zeros(3 * length(EBCN), 1);
    BCN3(1:3:end, 1) = 3 * EBCN - 2;
    BCN3(2:3:end, 1) = 3 * EBCN - 1;
    BCN3(3:3:end, 1) = 3 * EBCN - 0;

    element3DDofspzt = zeros(length(elementList), 3 * nElement3DNodes);

    element3DDofspzt(:, 1:3:end) = 3 * nodespzt(:, :) - 2;
    element3DDofspzt(:, 2:3:end) = 3 * nodespzt(:, :) - 1;
    element3DDofspzt(:, 3:3:end) = 3 * nodespzt(:, :) - 0;

    cpzt = 0;
    for elementNo = elementList % local matrices for each element in pzt
        cpzt = cpzt + 1;
        nodeList = SemMesh.element3DNodes(elementNo, :);

        [kufi, kfifi] = pzt_coupling_matrices(Materials.pztCoupling(:, :, iPiezoMat), ...
                                              Materials.pztPermittivity(:, :, iPiezoMat), SemMesh.nodeCoordinates(nodeList, 1), ...
                                              SemMesh.nodeCoordinates(nodeList, 2), SemMesh.nodeCoordinates(nodeList, 3), ...
                                              inverseVandermondeMatrixKsi, inverseVandermondeMatrixEta, ...
                                              inverseVandermondeMatrixDzeta, ksi, eta, dzeta, wKsi, wEta, wDzeta);
        kfiu = sparse(kufi');
        kfifi = sparse(kfifi);

        [iout] = dofs3d(element3DDofspzt, cpzt, SemMesh.nDofs, nElement3DNodes);
        [iout2] = dofs3dfifi(nodespzt, cpzt, SemMesh.nNodes, nElement3DNodes);

        Kfiu = Kfiu + iout2' * kfiu * iout; % assembly
        Kfifi = Kfifi + iout2' * kfifi * iout2; % assembly
    end
    PztSensor(iSensor).EBCI = EBCI;
    PztSensor(iSensor).EBCIN = [EBCI; EBCN];
    PztSensor(iSensor).EBCN = EBCN;
    PztSensor(iSensor).EBC0 = EBC0;
    PztSensor(iSensor).BC03 = BC03;
    PztSensor(iSensor).BCI3 = BCI3;
    PztSensor(iSensor).BCN3 = BCN3;
end

% selection of nodes and degress of freedom indexing for electric coupling matrices (all sensors)
sensorsEBCI = []; % internal nodes
sensorsEBCN = []; % top electrode nodes
% sensorsEBCIN = []; % internal nodes + top electrode nodes
sensorsBC03 = []; % grounded electrode degrees of freedom
sensorsBCI3 = []; % internal element degrees of freedom
sensorsBCN3 = []; % top electrode degrees of freedom

for iSensor = 1:Output.nSensors
    sensorsEBCI = [sensorsEBCI; PztSensor(iSensor).EBCI];
    sensorsEBCN = [sensorsEBCN; PztSensor(iSensor).EBCN];
    sensorsBC03 = [sensorsBC03; PztSensor(iSensor).BC03];
    sensorsBCI3 = [sensorsBCI3; PztSensor(iSensor).BCI3];
    sensorsBCN3 = [sensorsBCN3; PztSensor(iSensor).BCN3];

end
sensorsEBCIN = [sensorsEBCI; sensorsEBCN];
sensorsBC = [sensorsBC03; sensorsBCI3; sensorsBCN3];
KfifiO = [Kfifi(sensorsEBCI, sensorsEBCI), Kfifi(sensorsEBCI, sensorsEBCN)
          Kfifi(sensorsEBCN, sensorsEBCI), Kfifi(sensorsEBCN, sensorsEBCN)]; % sensor in open circuit
KfiuO = [Kfiu(sensorsEBCI, sensorsBC03), Kfiu(sensorsEBCI, sensorsBCI3), Kfiu(sensorsEBCI, sensorsBCN3)
         Kfiu(sensorsEBCN, sensorsBC03), Kfiu(sensorsEBCN, sensorsBCI3), Kfiu(sensorsEBCN, sensorsBCN3)];
invKfifiO = inv(KfifiO);
invKfifiO = full(invKfifiO);
if isGPUavailable  % push to gpu
    invKfifiO = gpuArray(invKfifiO);
    KfiuO = gpuArray(KfiuO);
end

clear element3DDofspzt Kfifi KfifiA Kfifin Kfiu;
% clear FiAi FiA;

%% Output files for solution
% extract response points from tags
%
Output.responseNodes = [];
counter = 0;
for i = 1:length(Output.responsePointTagList)
    for j = 1:length(SemMesh.PhysicalPoints)
        if Output.responsePointTagList(i) == SemMesh.PhysicalPoints(j).physicalTag
            counter = counter + 1;
            Output.responseNodes = [Output.responseNodes; SemMesh.PhysicalPoints(j).nodes];
            Output.responsePointPhysicalNames{counter} = SemMesh.PhysicalPoints(j).name;
        end
    end
end
Output.nResponseNodes = length(Output.responseNodes);

Output.outfileVoltage = fullfile(model_output_path, 'voltage');
Output.outfileDisplacements = fullfile(model_output_path, 'displacements');
Output.outfileVelocities = fullfile(model_output_path, 'velocities');
Output.outfileTime = fullfile(model_output_path, 'timeVector');

if Output.isVtuOutput
    [ParaviewMesh] = sem2paraview_mesh(SemMesh);
end
% save mesh related fields to one frame
if Output.isVtuOutput
    iSample = Output.sampleInterval;
    try
        save_vtu_frame_mesh(model_output_path, iSample, ParaviewMesh, Output.isBinary);
    catch
        if isVerbose
            if(~isempty(calculationProgress))
                calculationProgress.Value = 0;
                calculationProgress.Message = 'Problem saving mesh to vtu file';
                disp('Problem saving mesh to vtu file');
            else
                disp('Problem saving mesh to vtu file');
            end
        end
        isCalculationFailed = true;
        return
    end
    results_path =  fullfile(model_output_path, filesep, 'frames_vtu', filesep);
    file_name = ['frame', num2str(iSample, '%07u'), '.vtu'];
    source_file = fullfile(results_path, file_name);
    % copy the content to all remaining frames
    for iSample = 2 * Output.sampleInterval:Output.sampleInterval:Output.nFrames * Output.sampleInterval
        file_name = ['frame', num2str(iSample, '%07u'), '.vtu'];
        destination_file = fullfile(results_path, file_name);
        copyfile(source_file, destination_file);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U - displacements
% F - forces
% V - velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Initialize vectors';
    else
        disp('Initialize vectors');
    end
end
if isGPUavailable  % allocate on gpu
    Uold   = zeros(SemMesh.nDofs, 1, 'double', 'gpuArray');
    U      = zeros(SemMesh.nDofs, 1, 'double', 'gpuArray');
    Unew   = zeros(SemMesh.nDofs, 1, 'double', 'gpuArray');
    FXYZ_P = zeros(3 * nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    UX_P   = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    UY_P   = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    UZ_P   = zeros(nElements * nElement3DNodes, 1, 'double', 'gpuArray');
    V      = zeros(SemMesh.nDofs, 1, 'double', 'gpuArray');
else   % allocate on cpu
    Uold   = zeros(SemMesh.nDofs, 1);
    U      = zeros(SemMesh.nDofs, 1);
    Unew   = zeros(SemMesh.nDofs, 1);
    FXYZ_P = zeros(3 * nElements * nElement3DNodes, 1);
    UX_P   = zeros(nElements * nElement3DNodes, 1);
    UY_P   = zeros(nElements * nElement3DNodes, 1);
    UZ_P   = zeros(nElements * nElement3DNodes, 1);
    V      = zeros(SemMesh.nDofs, 1);
end

if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Initialize vectors. Done';
    else
        disp('Initialize vectors. Done');
    end
end
% Boundary conditions ??
% TODO: implement boundary conditions

%% Initial calculation - central difference scheme
a0 = 1 / (Excitation.dt^2);
a1 = 1 / (2 * Excitation.dt);
a2 = 2 * a0;
a3 = 1 / a2;
if isGPUavailable  % push to gpu
    a0 = gpuArray(a0);
    a1 = gpuArray(a1);
    a2 = gpuArray(a2);
    a3 = gpuArray(a3);
end
invEffMatrix = 1 ./ (a0 .* massMatrixGlobal + a1 * dampingMatrixGlobal);

if isGPUavailable  % allocate on gpu
    voltage = zeros(Excitation.nSamples, Output.nSensors, 'double', 'gpuArray');
    if Output.isDisplacementResponsePointsSelected
        displacementsX = zeros(Excitation.nSamples, Output.nResponseNodes, 'double', 'gpuArray');
        displacementsY = zeros(Excitation.nSamples, Output.nResponseNodes, 'double', 'gpuArray');
        displacementsZ = zeros(Excitation.nSamples, Output.nResponseNodes, 'double', 'gpuArray');
    end
    if Output.isVelocityResponsePointsSelected
        velocitiesX = zeros(Excitation.nSamples, Output.nResponseNodes, 'double', 'gpuArray');
        velocitiesY = zeros(Excitation.nSamples, Output.nResponseNodes, 'double', 'gpuArray');
        velocitiesZ = zeros(Excitation.nSamples, Output.nResponseNodes, 'double', 'gpuArray');
    end
else
    voltage = zeros(Excitation.nSamples, Output.nSensors);
    if Output.isDisplacementResponsePointsSelected
        displacementsX = zeros(Excitation.nSamples, Output.nResponseNodes);
        displacementsY = zeros(Excitation.nSamples, Output.nResponseNodes);
        displacementsZ = zeros(Excitation.nSamples, Output.nResponseNodes);
    end
    if Output.isVelocityResponsePointsSelected
        velocitiesX = zeros(Excitation.nSamples, Output.nResponseNodes);
        velocitiesY = zeros(Excitation.nSamples, Output.nResponseNodes);
        velocitiesZ = zeros(Excitation.nSamples, Output.nResponseNodes);
    end
end

timeFrames = zeros(Output.nFrames, 1);
timeElapsedSaveFrame = zeros(Output.nFrames, 1);
timeElapsedStep = zeros(Excitation.nSamples, 1);

if isGPUavailable  % allocate on gpu
    FiS = zeros(SemMesh.nNodes, 1, 'double', 'gpuArray'); % global vector of electric potential
else
    FiS = zeros(SemMesh.nNodes, 1); % global vector of electric potential
end
if isGPUavailable  % allocate on gpu
    KU_0 = zeros(SemMesh.nDofs, 1, 'double', 'gpuArray');
else  % allocate on cpu
    KU_0 = zeros(SemMesh.nDofs, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isVerbose
    if(~isempty(calculationProgress))
        calculationProgress.Value = 0;
        calculationProgress.Message = 'Time integration';
    else
        disp('Time integration');
    end
end
% profile on
c = 0;
for iSample = 2:Excitation.nSamples
    tstart = tic;
    % [iSample,Excitation.nSamples]
    % Check for Cancel button press
    if(~isempty(calculationProgress))
        if calculationProgress.CancelRequested
            isCalculationCanceled = true;
            return
        end
    end


    % transform from global displacements to local displacements stored by disconnected elements
    UX_P = U(indexMapX);
    UY_P = U(indexMapY);
    UZ_P = U(indexMapZ);

    % calculate strains
    [Exx, Eyy, Ezz, Exy, Eyz, Exz] = strains_spec_p(Nksi_P, Neta_P, Ndzeta_P, invJ11_P, invJ12_P, invJ13_P, invJ21_P, invJ22_P, invJ23_P, invJ31_P, invJ32_P, invJ33_P, UX_P, UY_P, UZ_P);

    % calculatate stresses
    [Sxx, Syy, Szz, Sxy, Syz, Sxz] = stresses_spec_p(Exx, Eyy, Ezz, Exy, Eyz, Exz, D11_P, D12_P, D13_P, D14_P, D22_P, D23_P, D24_P, D33_P, D34_P, D44_P, D55_P, D56_P, D66_P);

    % calculatate forces
    FX_P = NksiT_P * (Sxx .* invJ11_P .* WWWDetJ_P) + NetaT_P * (Sxx .* invJ21_P .* WWWDetJ_P) + NdzetaT_P * (Sxx .* invJ31_P .* WWWDetJ_P);
    FX_P = FX_P + NksiT_P * (Sxy .* invJ12_P .* WWWDetJ_P) + NetaT_P * (Sxy .* invJ22_P .* WWWDetJ_P) + NdzetaT_P * (Sxy .* invJ32_P .* WWWDetJ_P);
    FX_P = FX_P + NksiT_P * (Sxz .* invJ13_P .* WWWDetJ_P) + NetaT_P * (Sxz .* invJ23_P .* WWWDetJ_P) + NdzetaT_P * (Sxz .* invJ33_P .* WWWDetJ_P);

    FY_P = NksiT_P * (Syy .* invJ12_P .* WWWDetJ_P) + NetaT_P * (Syy .* invJ22_P .* WWWDetJ_P) + NdzetaT_P * (Syy .* invJ32_P .* WWWDetJ_P);
    FY_P = FY_P + NksiT_P * (Sxy .* invJ11_P .* WWWDetJ_P) + NetaT_P * (Sxy .* invJ21_P .* WWWDetJ_P) + NdzetaT_P * (Sxy .* invJ31_P .* WWWDetJ_P);
    FY_P = FY_P + NksiT_P * (Syz .* invJ13_P .* WWWDetJ_P) + NetaT_P * (Syz .* invJ23_P .* WWWDetJ_P) + NdzetaT_P * (Syz .* invJ33_P .* WWWDetJ_P);

    FZ_P = NksiT_P * (Szz .* invJ13_P .* WWWDetJ_P) + NetaT_P * (Szz .* invJ23_P .* WWWDetJ_P) + NdzetaT_P * (Szz .* invJ33_P .* WWWDetJ_P);
    FZ_P = FZ_P + NksiT_P * (Syz .* invJ12_P .* WWWDetJ_P) + NetaT_P * (Syz .* invJ22_P .* WWWDetJ_P) + NdzetaT_P * (Syz .* invJ32_P .* WWWDetJ_P);
    FZ_P = FZ_P + NksiT_P * (Sxz .* invJ11_P .* WWWDetJ_P) + NetaT_P * (Sxz .* invJ21_P .* WWWDetJ_P) + NdzetaT_P * (Sxz .* invJ31_P .* WWWDetJ_P);

    F = KU_0;  % clear global vector of internal forces

    FXYZ_P(1:3:end, :) = FX_P;
    FXYZ_P(2:3:end, :) = FY_P;
    FXYZ_P(3:3:end, :) = FZ_P;
    clear FX_P FY_P FZ_P;
    % assembly force vector - back to global vector
    for ibasket = 1:SemMesh.nBaskets
        F(I_G_dofs(:, ibasket)) = F(I_G_dofs(:, ibasket)) + FXYZ_P(I_L_dofs(:, ibasket));
    end
    if isGPUavailable  % allocate on gpu
        externalForces = zeros(SemMesh.nDofs, 1, 'double', 'gpuArray');
        topElectrodeVoltage = zeros(SemMesh.nNodes, 1, 'double', 'gpuArray');
        FiA = zeros(SemMesh.nNodes, 1, 'double', 'gpuArray');
    else
        externalForces = zeros(SemMesh.nDofs, 1);
        topElectrodeVoltage = zeros(SemMesh.nNodes, 1);
        FiA = zeros(SemMesh.nNodes, 1);
    end
    for iActuator = 1:Excitation.nActuators
        if Excitation.assignment(iActuator).category == 1  % Hann signal
            iHannSignal = Excitation.assignment(iActuator).index;
            topElectrodeVoltage(PztActuator(iActuator).EBCN, 1) = ...
                Excitation.Hann(iHannSignal).peakVoltage * Excitation.signals(iSample, iActuator);
        end
    end
    % internal charge
    fiAi = -invKfifiAi * (KfiuAi * U(actuatorsBC) +  KfifiAin * topElectrodeVoltage(actuatorsEBCN));
    % fiAi = - invKfifiAi*(KfifiAin * topElectrodeVoltage(actuatorsEBCN));

    FiA(actuatorsEBCN) = topElectrodeVoltage(actuatorsEBCN);
    FiA(actuatorsEBCI) = fiAi;
    FiA(actuatorsEBC0) = 0; % grounded electrodes

    % output signals
    % slowest line
    fio = -invKfifiO * (KfiuO * U(sensorsBC)); % electric potential in open circuit for all sensors
    FiS(sensorsEBCIN) = fio;
    % TODO: implement computation of electric field from the formula E = -Bfi*Fi
    for iSensor = 1:Output.nSensors
        voltage(iSample, iSensor) = -mean(FiS(PztSensor(iSensor).EBCN));
    end

    externalForces(actuatorsBC) = KfiuA' * FiA(actuatorsEBC); % local forces in pzt actuators

    % time integration
    F = -F + externalForces + (a2 .* massMatrixGlobal) .* U - a0 .* massMatrixGlobal .* Uold + a1 * dampingMatrixGlobal .* Uold; % piezoelectric excitation with damping

    Unew = invEffMatrix .* F;
    V = a1 * (-Uold + Unew);  % update velocity
    Uold = U;
    U = Unew;

    % additional pointwise outputs (displacements)
    if Output.isDisplacementResponsePointsSelected
        displacementsX(iSample, :) = U(3 * Output.responseNodes - 2, 1);
        displacementsY(iSample, :) = U(3 * Output.responseNodes - 1, 1);
        displacementsZ(iSample, :) =  U(3 * Output.responseNodes - 0, 1);
    end
    if Output.isVelocityResponsePointsSelected
        velocitiesX(iSample, :) = V(3 * Output.responseNodes - 2, 1);
        velocitiesY(iSample, :) = V(3 * Output.responseNodes - 1, 1);
        velocitiesZ(iSample, :) = V(3 * Output.responseNodes - 0, 1);
    end
    timeElapsedStep(iSample) = toc(tstart); % time step duration
    %% output frame
    % save frame to file
    if mod(iSample, Output.sampleInterval) == 0
        c = c + 1;
        timeFrames(c) = Excitation.timeVector(iSample);
        if max(abs(U)) > 100 * characteristicDistance || isnan(max(abs(U))) % any(isnan(U))
            if isVerbose
                if(~isempty(calculationProgress))
                    calculationProgress.Value = iSample / Excitation.nSamples;
                    calculationProgress.Message = 'Integration error';
                    disp('Integration error');
                else
                    disp('Integration error');
                end
                isCalculationFailed = true;
            end

            if isGPUavailable
                voltage = gather(voltage);
                if Output.isDisplacementResponsePointsSelected
                    displacementsX = gather(displacementsX);
                    displacementsY = gather(displacementsY);
                    displacementsZ = gather(displacementsZ);
                end
                if Output.isVelocityResponsePointsSelected
                    velocitiesX = gather(velocitiesX);
                    velocitiesY = gather(velocitiesY);
                    velocitiesZ = gather(velocitiesZ);
                end
            end
            save(Output.outfileVoltage, 'voltage');
            if Output.isDisplacementResponsePointsSelected
                save(Output.outfileDisplacements, 'displacementsX', 'displacementsY', 'displacementsZ');
            end
            if Output.isVelocityResponsePointsSelected
                save(Output.outfileVelocities, 'velocitiesX', 'velocitiesY', 'velocitiesZ');
            end
            timeVector = Excitation.timeVector;
            save(Output.outfileTime, 'timeVector', 'timeFrames');
            return
        end
        if Output.isDisplacementFramesSelected
            Uc = gather(U);
            frames_path = fullfile(model_output_path, filesep, 'frames_mat', filesep);
            if ~exist(frames_path, 'dir')
                mkdir(frames_path);
            end
            outfile_Ux = fullfile(model_output_path, filesep, 'frames_mat', filesep, ...
                                  ['Ux_frame', num2str(iSample)]);
            outfile_Uy = fullfile(model_output_path, filesep, 'frames_mat', filesep, ...
                                  ['Uy_frame', num2str(iSample)]);
            outfile_Uz = fullfile(model_output_path, filesep, 'frames_mat', filesep, ...
                                  ['Uz_frame', num2str(iSample)]);
            Ux = Uc(1:3:end);
            Uy = Uc(2:3:end);
            Uz = Uc(3:3:end);
            if Output.isMatOutput
                save(outfile_Ux, 'Ux', '-v7.3'); % frame output for global vector of Ux displacements
                save(outfile_Uy, 'Uy', '-v7.3'); % frame output for global vector of Uy displacements
                save(outfile_Uz, 'Uz', '-v7.3'); % frame output for global vector of Uz displacements
            end
        else
            Ux = [];
            Uy = [];
            Uz = [];
        end
        if Output.isVelocityFramesSelected
            Vc = gather(V);
            frames_path = fullfile(model_output_path, filesep, 'frames_mat', filesep);
            if ~exist(frames_path, 'dir')
                mkdir(frames_path);
            end
            outfile_Vx = fullfile(model_output_path, filesep, 'frames_mat', filesep, ...
                                  ['Vx_frame', num2str(iSample)]);
            outfile_Vy = fullfile(model_output_path, filesep, 'frames_mat', filesep, ...
                                  ['Vy_frame', num2str(iSample)]);
            outfile_Vz = fullfile(model_output_path, filesep, 'frames_mat', filesep, ...
                                  ['Vz_frame', num2str(iSample)]);
            Vx = Vc(1:3:end);
            Vy = Vc(2:3:end);
            Vz = Vc(3:3:end);
            if Output.isMatOutput
                save(outfile_Vx, 'Vx', '-v7.3'); % frame output for global vector of Vx velocities
                save(outfile_Vy, 'Vy', '-v7.3'); % frame output for global vector of Vy velocities
                save(outfile_Vz, 'Vz', '-v7.3'); % frame output for global vector of Vz velocities
            end
        else
            Vx = [];
            Vy = [];
            Vz = [];
        end
        if Output.isVtuOutput
            if Output.isDisplacementFramesSelected || Output.isVelocityFramesSelected
                try
                    save_vtu_frame_data(model_output_path, iSample, ...
                                   Ux, Uy, Uz, Vx, Vy, Vz, Output.isBinary);
                catch
                    if isVerbose
                        disp('Problem saving data to vtu file');
                    end
                    isCalculationFailed = true;
                    return
                end
            end
        end
        if(~isempty(SignalPreview))
            plot(SignalPreview, Excitation.timeVector * 1e3, voltage(:, Output.sensorPreviewNo));
        end
        drawnow;
        % get time elapsed for frame saving
        timeElapsedSaveFrame(c) = toc(tstart) - timeElapsedStep(iSample);   
        averageSaveFrameTime = mean(timeElapsedSaveFrame(1:c));
        averageTimeStep = mean(timeElapsedStep(2:iSample));
        timeElapsed = averageSaveFrameTime + averageTimeStep * Output.sampleInterval;   
        expectedDurationSeconds = seconds((Output.nFrames - c) * timeElapsed);
        progress = round(iSample / Excitation.nSamples * 100);

        timeNow = datetime('now', 'TimeZone', 'local', 'Format', 'd-MMM-y HH:mm');
        timeExpected = timeNow + expectedDurationSeconds;
        message = sprintf('Time integration\nTask progress: %d%%\nExpected finish time:\n%s', ...
                          progress, string(timeExpected));
        if isVerbose
            if(~isempty(calculationProgress))
                calculationProgress.Value = iSample / Excitation.nSamples;
                calculationProgress.Message = message;
            else
                disp(message);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  END OF MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Output.isVtuOutput
    save_pvd_animation(model_output_path, Output.nFrames, Output.sampleInterval, Excitation.timeVector);
end

if isGPUavailable
    voltage = gather(voltage);
    if Output.isDisplacementResponsePointsSelected
        displacementsX = gather(displacementsX);
        displacementsY = gather(displacementsY);
        displacementsZ = gather(displacementsZ);
    end
    if Output.isVelocityResponsePointsSelected
        velocitiesX = gather(velocitiesX);
        velocitiesY = gather(velocitiesY);
        velocitiesZ = gather(velocitiesZ);
    end
end
save(Output.outfileVoltage, 'voltage');
if Output.isDisplacementResponsePointsSelected
    save(Output.outfileDisplacements, 'displacementsX', 'displacementsY', 'displacementsZ');
end
if Output.isVelocityResponsePointsSelected
    save(Output.outfileVelocities, 'velocitiesX', 'velocitiesY', 'velocitiesZ');
end
timeVector = Excitation.timeVector;
save(Output.outfileTime, 'timeVector', 'timeFrames');

% ---------------------------------------------------------------------------------------------------

end
