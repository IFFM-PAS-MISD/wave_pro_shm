function [Exx, Eyy, Ezz, gammaxy, gammayz, gammaxz] = strains_spec_p(Nksi,   Neta, Ndzeta, ...
                                                                   invJ11, invJ12, invJ13, ...
                                                                   invJ21, invJ22, invJ23, ...
                                                                   invJ31, invJ32, invJ33, ...
                                                                       UX,     UY,     UZ)
%
% calculate strains at nodes of spectral elements
% suitable for parallel mode
%
% USAGE::
%
%   [Exx,Eyy,Ezz,gammaxy,gammayz,gammaxz] = strains_spec_p(Nksi,Neta,Ndzeta,...
%                                              invJ11,invJ12,invJ13,...
%                                              invJ21,invJ22,invJ23, ...
%                                              invJ31,invJ32,invJ33, ...
%                                              UX,UY,UZ)
%
% Arguments:
%     Nksi (double):
%       sparse block-diagonal matrix of shape function derivatives along ksi axis
%
%     Neta (double):
%       sparse block-diagonal matrix of shape function derivatives along ksi axis
%
%     Ndzeta (double):
%       sparse block-diagonal matrix of shape function derivatives along dzeta axis
%
%
%     invJij (double):
%       inverse Jacobian of ij element, i,j = 1,2,3 (vectors)
%
%     Ux, Uy, Uz (double):
%       global displacements at element nodes; dimensions are the same
%       as for inverse Jacobians invJij
%
%
% Returns:
%     Exx, Eyy, Ezz, gammaxy, gammayz, gammaxz (double):
%       strains at element nodes (vectors)
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

% note - brackets are important!

Exx = (Nksi * UX) .* invJ11 + (Neta * UX) .* invJ21 + (Ndzeta * UX) .* invJ31;

Eyy = (Nksi * UY) .* invJ12 + (Neta * UY) .* invJ22 + (Ndzeta * UY) .* invJ32;

Ezz = (Nksi * UZ) .* invJ13 + (Neta * UZ) .* invJ23 + (Ndzeta * UZ) .* invJ33;

B41UX = (Nksi * UX) .* invJ12 + (Neta * UX) .* invJ22 + (Ndzeta * UX) .* invJ32;
B42UY = (Nksi * UY) .* invJ11 + (Neta * UY) .* invJ21 + (Ndzeta * UY) .* invJ31;

gammaxy = B41UX + B42UY;

B52UY = (Nksi * UY) .* invJ13 + (Neta * UY) .* invJ23 + (Ndzeta * UY) .* invJ33;
B53UZ = (Nksi * UZ) .* invJ12 + (Neta * UZ) .* invJ22 + (Ndzeta * UZ) .* invJ32;

gammayz = B52UY + B53UZ;

B61UX = (Nksi * UX) .* invJ13 + (Neta * UX) .* invJ23 + (Ndzeta * UX) .* invJ33;
B63UZ = (Nksi * UZ) .* invJ11 + (Neta * UZ) .* invJ21 + (Ndzeta * UZ) .* invJ31;

gammaxz = B61UX + B63UZ;

% ---------------------------------------------------------------------------------------------------

end
