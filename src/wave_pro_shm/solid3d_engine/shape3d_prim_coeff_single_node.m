function [Nx, Ny, Nz] = shape3d_prim_coeff_single_node(nx, ny, nz, ...
                                                     inverseVandermondeMatrixKsi, ...
                                                     inverseVandermondeMatrixEta, ...
                                                     inverseVandermondeMatrixDzeta, ksi, eta, dzeta)
%
% Calculate shape function derivatives at single Gauss-Lobatto-Legendre node
% of a 3D solid spectral element
%
% USAGE::
%
%   [Nx,Ny,Nz] = shape3d_prim_coeff_single_node(nx,ny,nz, ...
%                                               inverseVandermondeMatrixKsi,...
%                                               inverseVandermondeMatrixEta, ...
%                                               inverseVandermondeMatrixDzeta,ksi,eta,dzeta)
%
% Arguments:
%     nx, ny, nz (integer):
%       number of nodes along ksi, eta, and dzeta, respectively
%
%     inverseVandermondeMatrixKsi (double):
%       inverse Vandermonde Matrix along ksi axis
%
%     inverseVandermondeMatrixEta (double):
%       inverse Vandermonde Matrix along eta axis
%
%     inverseVandermondeMatrixDzeta (double):
%       inverse Vandermonde Matrix along dzeta axis
%
%     ksi,eta,dzeta (double):
%       local coordinates at selected node of 3D solid spectral element
%
% Returns:
%     Nx, Ny, Nz (double):
%       shape function derivatives in respect to ksi, eta and dzeta
%       dimensions [1,nx*ny*nz]
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

Nx = zeros(1, nx * ny * nz);
Ny = zeros(1, nx * ny * nz);
Nz = zeros(1, nx * ny * nz);
ux = zeros(1, nx);
uy = zeros(1, ny);
uz = zeros(1, nz);

ux(1) = 1;
uy(1) = 1;
uz(1) = 1;

for j = 2:nx
    ux(j) = ksi^(j - 1);
end
for j = 2:ny
    uy(j) = eta^(j - 1);
end
for j = 2:nz
    uz(j) = dzeta^(j - 1);
end
Sx = ux * inverseVandermondeMatrixKsi';
Sy = uy * inverseVandermondeMatrixEta';
Sz = uz * inverseVandermondeMatrixDzeta';

uxprim = zeros(1, nx);
uyprim = zeros(1, ny);
uzprim = zeros(1, nz);

for j = 2:nx
    uxprim(j) = (j - 1) * ksi^(j - 2);
end
for j = 2:ny
    uyprim(j) = (j - 1) * eta^(j - 2);
end
for j = 2:nz
    uzprim(j) = (j - 1) * dzeta^(j - 2);
end
Sxprim = inverseVandermondeMatrixKsi * uxprim';
Syprim = uyprim * inverseVandermondeMatrixEta';
Szprim = uzprim * inverseVandermondeMatrixDzeta';
cc = 0;
for j3 = 1:nz
    for j2 = 1:ny % shape function number
        for j1 = 1:nx
            cc = cc + 1;
            Nx(cc) = Sxprim(j1) * Sy(j2) * Sz(j3);
            if abs(Nx(cc)) < 1e-12 Nx(cc) = 0;
            end
            Ny(cc) = Sx(j1) * Syprim(j2) * Sz(j3);
            if abs(Ny(cc)) < 1e-12 Ny(cc) = 0;
            end
            Nz(cc) = Sx(j1) * Sy(j2) * Szprim(j3);
            if abs(Nz(cc)) < 1e-12 Nz(cc) = 0;
            end
        end
    end
end
% ---------------------------------------------------------------------------------------------------

end
