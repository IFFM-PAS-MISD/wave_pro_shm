function [elasticityMatrixPiezo] = elasticity_matrix_piezo(YE11, YE33, YE55, nu12)
%
% Elasticity matrix of piezoelectric material (transversly isotropic)
%
% USAGE::
%
%   [elasticityMatrixPiezo] = elasticity_matrix_piezo(YE11,YE33,YE55,nu12)
%
% Arguments:
%     YE11 (double):
%       Young modulus in principal direction 11, Units [Pa]
%
%     YE33 (double):
%       Young modulus in principal direction 33, Units [Pa]
%
%     YE55 (double):
%       Young modulus in principal direction 55, Units [Pa]
%
%     nu12 (double):
%       Poisson ratio
%
% Returns:
%     elasticityMatrixPiezo (double):
%       Elasticity matrix, dimensions [6, 6], Units [Pa]
%
% .. Note:: Notation Sxx,Syy,Szz,Sxy,Syz,Sxz is different then Voigt notation!
%
%    Transverse isotropy ::
%
%    {Sxx}    |D11 D12 D13   0   0   0| {    Exx}
%    {Syy}    |D12 D22 D13   0   0   0| {    Eyy}
%    {Szz} =  |D13 D12 D33   0   0   0| {    Ezz}
%    {Sxy}    |  0   0   0 D66   0   0| {gammaxy}
%    {Syz}    |  0   0   0   0 D44   0| {gammayz}
%    {Sxz}    |  0   0   0   0   0 D55| {gammaxz}
%
% .. math::
%    \begin{Bmatrix}
%    \sigma_{xx} \\ \sigma_{yy} \\ \sigma_{zz} \\ \sigma_{xy} \\ \sigma_{yz} \\ \sigma_{xz} \end{Bmatrix} =
%    \begin{bmatrix}
%    D_{11} & D_{12} & D_{13} &      0 &      0 &      0\\
%    D_{12} & D_{22} & D_{13} &      0 &      0 &      0\\
%    D_{13} & D_{12} & D_{33} &      0 &      0 &      0\\
%         0 &      0 &      0 & D_{66} &      0 &      0\\
%         0 &      0 &      0 &      0 & D_{44} &      0\\
%         0 &      0 &      0 &      0 &      0 & D_{55}
%    \end{bmatrix}
%    \begin{Bmatrix}
%    \epsilon_{xx} \\ \epsilon_{yy} \\ \epsilon_{zz} \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{xz} \end{Bmatrix}
%
% .. seealso::
%
%    * https://support.onscale.com/hc/en-us/articles/360002073378-Calculating-Piezoelectric-Material-Properties-from-Material-Datasheet
%
%    * A. Deraemaeker, S. Benelechi, A. Benjeddou, and A. Preumont.
%      Analytical and Numerical Computation of Homogenized Properties of MFCS:
%      Application to Composite Boom with MFC Actuators and Sensors. In W. Ostachowicz,
%      J. Holnicki-szulc, and C. M. Soares, editors, *III ECCOMAS Thematic Conference
%      on Smart Structures and Materials*, 1–22. Gdansk, Poland, July 9-11, 2007.
%      European Community on Computational Methods in Applied Science - ECCOMAS.
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------
YE22 = YE11;
G12 = YE11 / (2 * (1 + nu12));

% elastic compliance matrix
Spzt = [  1 / YE11, -nu12 / YE22, -nu12 / YE33,       0,          0,         0;
      -nu12 / YE22,     1 / YE22, -nu12 / YE33,       0,          0,         0;
      -nu12 / YE33, -nu12 / YE33,     1 / YE33,       0,          0,         0;
                 0,            0,            0, 1 / G12,          0,         0;
                 0,            0,            0,       0, 1 / (YE55),         0;
                 0,            0,            0,       0,          0, 1 / (YE55)];

if YE11 == 0 || YE33 == 0 || YE55 == 0
    elasticityMatrixPiezo = zeros(6, 6);
else
    elasticityMatrixPiezo = inv(Spzt);
end
% ---------------------------------------------------------------------------------------------------

end
