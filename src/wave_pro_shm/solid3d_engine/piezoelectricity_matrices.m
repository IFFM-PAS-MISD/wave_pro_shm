function [chargeConstantMatrixd,permittivityMatrixEpsS,chargeConstantMatrixe] = ...
          piezoelectricity_matrices(chargeConstantd31, chargeConstantd33, chargeConstantd15, ...
          relativePermittivityEpsrS11, relativePermittivityEpsrS33, elasticityMatrixPZT)
    % 
    % Permittivity, charge and voltage constant matrices
    % 
    % USAGE:: 
    %
    %   [chargeConstantMatrixd,permittivityMatrixEpsS,chargeConstantMatrixe] = ...
    %    piezoelectricity_matrices(chargeConstantd31,chargeConstantd33,chargeConstantd15,...
    %    relativePermittivityEpsrS11,relativePermittivityEpsrS33,elasticityMatrixPZT)
    % 
    % Arguments: 
    %     chargeConstantd31 (double): 
    %       charge constant, Units [C/N] or [m/V] 
    % 
    %     chargeConstantd33 (double): 
    %       charge constant, Units [C/N] or [m/V] 
    % 
    %     chargeConstantd15 (double): 
    %       charge constant, Units [C/N] or [m/V] 
    %
    %     relativePermittivityEpsrS11 (double): 
    %       relative electric permittivity (dielectric constant) at zero strain
    %
    %     relativePermittivityEpsrS33 (double): 
    %       relative electric permittivity (dielectric constant) at zero strain
    %
    %     elasticityMatrixPZT (double): 
    %       elasticity Matrix, dimensions [6,6], Units [Pa]
    % 
    % Returns: 
    %     chargeConstantMatrixd (double): 
    %       charge constant matrix d, dimensions [3, 6], Units [C/N]
    %
    %     permittivityMatrixEpsS (double): 
    %       permittivity (dielectric) matrix EpsS, dimensions [3, 3], Units [Fa/m]
    %
    %     chargeConstantMatrixe (double): 
    %       charge constant matrix e, dimensions [3, 6], Units [C/m^2]
    %
    % .. Note:: Notation Sxx,Syy,Szz,Sxy,Syz,Sxz is diffrent than Voigt notation! 
    % 
    % 
    % (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl 
    % Institute of Fluid Flow Machinery, Polish Academy of Sciences 
    % Mechanics of Intelligent Structures Department 

%--------------------------------------------------------------------------------------------------- 

eps0 = 8.854e-12; % [F/m] vacuum permittivity
% electric permittivity (dielectric constant) at zero stress, Units: [F/m]
epsS11 = relativePermittivityEpsrS11 * eps0; 
epsS33 = relativePermittivityEpsrS33 * eps0; 
permittivityMatrixEpsS = [epsS11,      0,      0;
                               0, epsS11,      0;
                               0,      0, epsS33];

                                       %sx                 sy                 sz                sxy                sxz                syz         
chargeConstantMatrixd = [                0,                 0,                 0,                 0, chargeConstantd15,                 0;
                                         0,                 0,                 0,                 0,                 0, chargeConstantd15;
                         chargeConstantd31, chargeConstantd31, chargeConstantd33,                 0,                 0,                0];




chargeConstantMatrixe = chargeConstantMatrixd * elasticityMatrixPZT;
%---------------------------------------------------------------------------------------------------

end

