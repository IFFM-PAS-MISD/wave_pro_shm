function [elasticityMatrix] = elasticity_matrix_isotropy(youngModulus, poissonRatio)
%
% Elasticity matrix of isotropic material
%
% USAGE::
%
%   [elasticityMatrix] = elasticity_matrix_isotropy(youngModulus,poissonRatio)
%
% Arguments:
%     youngModulus (double):
%       Young modulus, Units [Pa]
%
%     poissonRatio (double):
%       Poisson ratio
%
% Returns:
%     elasticityMatrix (double):
%       Elasticity matrix, dimensions [6, 6], Units [Pa]
%
%
% .. Note:: Notation Sxx,Syy,Szz,Sxy,Syz,Sxz is different then Voigt notation! ::
%
%    {Sxx}    |D11 D12 D12   0   0   0| {    Exx}
%    {Syy}    |D12 D11 D12   0   0   0| {    Eyy}
%    {Szz} =  |D12 D12 D12   0   0   0| {    Ezz}
%    {Sxy}    |  0   0   0 D44   0   0| {gammaxy}
%    {Syz}    |  0   0   0   0 D44   0| {gammayz}
%    {Sxz}    |  0   0   0   0   0 D44| {gammaxz}
%
%
% .. math::
%
%    \begin{Bmatrix}
%    \sigma_{xx} \\ \sigma_{yy} \\ \sigma_{zz} \\ \sigma_{xy} \\ \sigma_{yz} \\ \sigma_{xz}
%    \end{Bmatrix} =
%    \begin{bmatrix}
%    D_{11} & D_{12} & D_{12} &      0 &      0 &      0\\
%    D_{12} & D_{11} & D_{12} &      0 &      0 &      0\\
%    D_{12} & D_{12} & D_{12} &      0 &      0 &      0\\
%         0 &      0 &      0 & D_{44} &      0 &      0\\
%         0 &      0 &      0 &      0 & D_{44} &      0\\
%         0 &      0 &      0 &      0 &      0 & D_{44}
%    \end{bmatrix}
%    \begin{Bmatrix}
%    \epsilon_{xx} \\ \epsilon_{yy} \\ \epsilon_{zz} \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{xz}
%    \end{Bmatrix}
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

stiffnessCoeff = youngModulus / ((1 + poissonRatio) * (1 - 2 * poissonRatio));

elasticityMatrix = [1 - poissonRatio,     poissonRatio,     poissonRatio,                          0,                          0,                          0;
                        poissonRatio, 1 - poissonRatio,     poissonRatio,                          0,                          0,                          0;
                        poissonRatio,     poissonRatio, 1 - poissonRatio,                          0,                          0,                          0;
                                   0,                0,                0, (1 - 2 * poissonRatio) / 2,                          0,                          0;
                                   0,                0,                0,                          0, (1 - 2 * poissonRatio) / 2,                          0;
                                   0,                0,                0,                          0,                          0, (1 - 2 * poissonRatio) / 2] * stiffnessCoeff;

% ---------------------------------------------------------------------------------------------------

end
