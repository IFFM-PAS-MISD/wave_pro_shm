function [Nx, Ny, Nz] = shape3d_prim_coeff(nx, ny, nz, inverseVandermondeMatrixKsi, ...
                                         inverseVandermondeMatrixEta, ...
                                         inverseVandermondeMatrixDzeta, ksi, eta, dzeta)
%
% Calculate shape function derivatives at Gauss-Lobatto-Legendre nodes of a 3D solid spectral
% element
%
% USAGE::
%
%   [Nx,Ny,Nz] = shape3d_prim_coeff(nx,ny,nz,inverseVandermondeMatrixKsi,...
%                                   inverseVandermondeMatrixEta, ...
%                                   inverseVandermondeMatrixDzeta,ksi,eta,dzeta)
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
%       local coordinates of nodes of 3D solid spectral element (vectors)
%       dimensions ksi [nx], eta [ny], dzeta [nz]
%
% Returns:
%     Nx, Ny, Nz (double):
%       shape function derivatives in respect to ksi, eta and dzeta
%       dimensions [nx*ny*nz,nx*ny*nz]
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

Nx = zeros(nx * ny * nz, nx * ny * nz);
Ny = zeros(nx * ny * nz, nx * ny * nz);
Nz = zeros(nx * ny * nz, nx * ny * nz);
ux = zeros(nx, nx);
uy = zeros(ny, ny);
uz = zeros(nz, nz);

ux(:, 1) = 1;
uy(:, 1) = 1;
uz(:, 1) = 1;

for j = 2:nx
    ux(:, j) = ksi.^(j - 1);
end
for j = 2:ny
    uy(:, j) = eta.^(j - 1);
end
for j = 2:nz
    uz(:, j) = dzeta.^(j - 1);
end
Sx = ux * inverseVandermondeMatrixKsi';
Sy = uy * inverseVandermondeMatrixEta';
Sz = uz * inverseVandermondeMatrixDzeta';

uxprim = zeros(nx, nx);
uyprim = zeros(ny, ny);
uzprim = zeros(nz, nz);

for j = 2:nx
    uxprim(:, j) = (j - 1) * ksi.^(j - 2);
end
for j = 2:ny
    uyprim(:, j) = (j - 1) * eta.^(j - 2);
end
for j = 2:nz
    uzprim(:, j) = (j - 1) * dzeta.^(j - 2);
end
Sxprim = uxprim * inverseVandermondeMatrixKsi';
Syprim = uyprim * inverseVandermondeMatrixEta';
Szprim = uzprim * inverseVandermondeMatrixDzeta';
c = 0;
for k3 = 1:nz
    for k2 = 1:ny
        for k1 = 1:nx
            cc = 0;
            c = c + 1;
            for j3 = 1:nz
                for j2 = 1:ny % shape function number
                    for j1 = 1:nx
                        cc = cc + 1;
                        Nx(c, cc) = Sxprim(k1, j1) * Sy(k2, j2) * Sz(k3, j3);
                        if abs(Nx(c, cc)) < 1e-12
                            Nx(c, cc) = 0;
                        end
                        Ny(c, cc) = Sx(k1, j1) * Syprim(k2, j2) * Sz(k3, j3);
                        if abs(Ny(c, cc)) < 1e-12
                            Ny(c, cc) = 0;
                        end
                        Nz(c, cc) = Sx(k1, j1) * Sy(k2, j2) * Szprim(k3, j3);
                        if abs(Nz(c, cc)) < 1e-12
                            Nz(c, cc) = 0;
                        end
                    end
                end
            end
        end
    end
end
% ---------------------------------------------------------------------------------------------------

end
