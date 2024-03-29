function [Sxx, Syy, Szz, Sxy, Syz, Sxz] = stresses_spec_p(Exx, Eyy, Ezz, gammaxy, gammayz, gammaxz, ...
                                                          D11, D12, D13, D14, D22, D23, D24, D33, ...
                                                          D34, D44, D55, D56, D66)
%
% calculate stresses at nodes of spectral elements
% suitable for parallel mode
%
% USAGE::
%
%   [Sxx,Syy,Szz,Sxy,Syz,Sxz] = stresses_spec_p(Exx,Eyy,Ezz,gammaxy,gammayz,gammaxz, ...
%                               D11,D12,D13,D14,D22,D23,D24,D33,D34,D44,D55,D56,D66)
%
% Arguments:
%     Exx, Eyy, Ezz, gammaxy, gammayz, gammaxz (double):
%       strains at element nodes (vectors)
%
%     Dij (double):
%       elasticity matrix coefficients (vectors); i,j = 1:6
%
% Returns:
%     Sxx, Syy, Szz, Sxy, Syz, Sxz (double):
%       stresses at element nodes (vectors)
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

Sxx = Exx .* D11 + Eyy .* D12 + Ezz .* D13 + gammaxy .* D14;
Syy = Exx .* D12 + Eyy .* D22 + Ezz .* D23 + gammaxy .* D24;
Szz = Exx .* D13 + Eyy .* D23 + Ezz .* D33 + gammaxy .* D34;
Sxy = Exx .* D14 + Eyy .* D24 + Ezz .* D34 + gammaxy .* D44;
Syz = gammayz .* D55 + gammaxz .* D56;
Sxz = gammayz .* D56 + gammaxz .* D66;

% ---------------------------------------------------------------------------------------------------

end
