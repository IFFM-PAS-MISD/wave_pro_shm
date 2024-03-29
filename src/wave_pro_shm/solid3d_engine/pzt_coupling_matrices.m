function [kufi, kfifi] = pzt_coupling_matrices(chargeConstantMatrixe, permittivityMatrixEpsS, x, y, z, ...
                                              inverseVandermondeMatrixKsi, ...
                                              inverseVandermondeMatrixEta, ...
                                              inverseVandermondeMatrixDzeta, ...
                                              ksi, eta, dzeta, wKsi, wEta, wDzeta)
%
% Calculate piezoelectric coupling matrices of a single 3D solid spectral element
%
% USAGE::
%
%   [kufi,kfifi] = pzt_coupling_matrices(chargeConstantMatrixe,permittivityMatrixEpsS,x,y,z,...
%                                        inverseVandermondeMatrixKsi,...
%                                        inverseVandermondeMatrixEta,...
%                                        inverseVandermondeMatrixDzeta,...
%                                        ksi,eta,dzeta,wKsi,wEta,wDzeta)
%
% Arguments:
%     chargeConstantMatrixe (double):
%       charge constant matrix e (stress-charge form), dimensions [3, 6], Units [C/m^2]
%
%     permittivityMatrixEpsS (double):
%       permittivity (dielectric) matrix EpsS, dimensions [3, 3], Units [Fa/m]
%
%     x,y,z (double):
%       coordinates of nodes of 3D solid spectral element
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
%     ksi, eta, dzeta (double):
%       local coordinates in range -1:1 (Gauss-Lobatto-Legendre points),
%       row vector of dimension dependent on the approximation order
%
%     wKsi,wEta,wDzeta (double):
%       weights at Gauss-Lobatto-Legendre points
%
% Returns:
%     kufi (double):
%       displacement-charge matrix coupling
%
%     kfifi (double):
%       charge matrix coupling
%
%
% .. seealso:: Function :func:`shape3d_prim_coeff_single_node`
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

nx = length(ksi);
ny = length(eta);
nz = length(dzeta);

% 3D case
B = zeros(6, nx * ny * nz * 3);
Bfi = zeros(3, nx * ny * nz);
kufi = zeros(nx * ny * nz * 3, nx * ny * nz);
kfifi = zeros(nx * ny * nz, nx * ny * nz);
for kz = 1:nz % dzeta
    cm = 0;
    for ky = 1:ny % eta
        for kx = 1:nx % ksi
            J = zeros(3, 3);
            kk = 0;
            [Nprimx, Nprimy, Nprimz] = shape3d_prim_coeff_single_node(nx, ny, nz, ...
                                                                    inverseVandermondeMatrixKsi, ...
                                                                    inverseVandermondeMatrixEta, ...
                                                                    inverseVandermondeMatrixDzeta, ...
                                                                    ksi(kx), eta(ky), dzeta(kz));
            for j3 = 1:nz
                for j2 = 1:ny % shape function number
                    for j1 = 1:nx
                        kk = kk + 1;
                        J(1, 1) = J(1, 1) + Nprimx(kk) * x(kk);
                        J(1, 2) = J(1, 2) + Nprimx(kk) * y(kk);
                        J(1, 3) = J(1, 3) + Nprimx(kk) * z(kk);

                        J(2, 1) = J(2, 1) + Nprimy(kk) * x(kk);
                        J(2, 2) = J(2, 2) + Nprimy(kk) * y(kk);
                        J(2, 3) = J(2, 3) + Nprimy(kk) * z(kk);

                        J(3, 1) = J(3, 1) + Nprimz(kk) * x(kk);
                        J(3, 2) = J(3, 2) + Nprimz(kk) * y(kk);
                        J(3, 3) = J(3, 3) + Nprimz(kk) * z(kk);
                    end
                end
            end
            %          if(abs(J(1,1)<eps)) J(1,1)=0; end;
            %          if(abs(J(1,2)<eps)) J(1,2)=0; end;
            %          if(abs(J(1,3)<eps)) J(1,3)=0; end;
            %          if(abs(J(2,1)<eps)) J(2,1)=0; end;
            %          if(abs(J(2,2)<eps)) J(2,2)=0; end;
            %          if(abs(J(2,3)<eps)) J(2,3)=0; end;
            %          if(abs(J(3,1)<eps)) J(3,1)=0; end;
            %          if(abs(J(3,2)<eps)) J(3,2)=0; end;
            %          if(abs(J(3,3)<eps)) J(3,3)=0; end;
            Jj = [J(1, 1), J(1, 2), J(1, 3), J(2, 1), J(2, 2), J(2, 3), J(3, 1), J(3, 2), J(3, 3)];

            % normal vector
            V3(1) = Jj(2) * Jj(6) - Jj(5) * Jj(3);
            V3(2) = Jj(4) * Jj(3) - Jj(1) * Jj(6);
            V3(3) = Jj(1) * Jj(5) - Jj(4) * Jj(2);
            y0 = V3(1) * V3(1) + V3(2) * V3(2) + V3(3) * V3(3);
            ys = sqrt(y0);
            v31 = V3(1) / ys;
            v32 = V3(2) / ys;
            v33 = V3(3) / ys;
            % /*v1 vector*/
            V1(1) = Jj(1);
            V1(2) = Jj(2);
            V1(3) = Jj(3);
            y0 = V1(1) * V1(1) + V1(2) * V1(2) + V1(3) * V1(3);
            ys = sqrt(y0);
            v11 = V1(1) / ys;
            v12 = V1(2) / ys;
            v13 = V1(3) / ys;
            % /*cross product v2=v3xv1*/
            v21 = v32 * v13 - v33 * v12;
            v22 = -v31 * v13 + v33 * v11;
            v23 = v31 * v12 - v32 * v11;
            y0 = v21 * v21 + v22 * v22 + v23 * v23;
            ys = sqrt(y0);
            v21 = v21 / ys;
            v22 = v22 / ys;
            v23 = v23 / ys;
            %%%
            % /* determinant and inverse Jacobian  */
            inJ = inv(J);
            detJ = det(J);
            wwwdetJ = wKsi(kx) * wEta(ky) * wDzeta(kz) * detJ;
            cc = 0;
            for j3 = 1:nz
                for j2 = 1:ny % shape function number
                    for j1 = 1:nx
                        cc = cc + 1;
                        Nx = Nprimx(cc) * inJ(1, 1) + Nprimy(cc) * inJ(1, 2) + Nprimz(cc) * inJ(1, 3);
                        Ny = Nprimx(cc) * inJ(2, 1) + Nprimy(cc) * inJ(2, 2) + Nprimz(cc) * inJ(2, 3);
                        Nz = Nprimx(cc) * inJ(3, 1) + Nprimy(cc) * inJ(3, 2) + Nprimz(cc) * inJ(3, 3);
                        % B - strain - nodal displacement matrix;
                        B11 = Nx * v11^2 + Ny * v11 * v12 + Nz * v11 * v13;
                        B12 = Nx * v11 * v12 + Ny * v12^2 + Nz * v12 * v13;
                        B13 = Nx * v11 * v13 + Ny * v12 * v13 + Nz * v13^2;
                        B21 = Nx * v21^2 + Ny * v21 * v22 + Nz * v21 * v23;
                        B22 = Nx * v21 * v22 + Ny * v22^2 + Nz * v22 * v23;
                        B23 = Nx * v21 * v23 + Ny * v22 * v23 + Nz * v23^2;
                        B31 = Nx * v31^2 + Ny * v31 * v32 + Nz * v31 * v33;
                        B32 = Nx * v31 * v32 + Ny * v32^2 + Nz * v32 * v33;
                        B33 = Nx * v31 * v33 + Ny * v32 * v33 + Nz * v33^2;
                        B41 = 2 * Nx * v11 * v21 + Ny * (v12 * v21 + v11 * v22) + Nz * (v13 * v21 + v11 * v23);
                        B42 = 2 * Ny * v12 * v22 + Nx * (v12 * v21 + v11 * v22) + Nz * (v13 * v22 + v12 * v23);
                        B43 = 2 * Nz * v13 * v23 + Nx * (v13 * v21 + v11 * v23) + Ny * (v13 * v22 + v12 * v23);
                        B51 = 2 * Nx * v21 * v31 + Ny * (v22 * v31 + v21 * v32) + Nz * (v23 * v31 + v21 * v33);
                        B52 = 2 * Ny * v22 * v32 + Nx * (v22 * v31 + v21 * v32) + Nz * (v23 * v32 + v22 * v33);
                        B53 = 2 * Nz * v23 * v33 + Nx * (v23 * v31 + v21 * v33) + Ny * (v23 * v32 + v22 * v33);
                        B61 = 2 * Nx * v11 * v31 + Ny * (v12 * v31 + v11 * v32) + Nz * (v13 * v31 + v11 * v33);
                        B62 = 2 * Ny * v12 * v32 + Nx * (v12 * v31 + v11 * v32) + Nz * (v13 * v32 + v12 * v33);
                        B63 = 2 * Nz * v13 * v33 + Nx * (v13 * v31 + v11 * v33) + Ny * (v13 * v32 + v12 * v33);

                        % epsx=du/dx
                        B(1, 3 * (cc) - 2) = B11;
                        B(1, 3 * (cc) - 1) = B12;
                        B(1, 3 * (cc)) = B13;
                        % epsy=dv/dy
                        B(2, 3 * (cc) - 2) = B21;
                        B(2, 3 * (cc) - 1) = B22;
                        B(2, 3 * (cc)) = B23;
                        % epsz=dw/dz
                        B(3, 3 * (cc) - 2) = B31;
                        B(3, 3 * (cc) - 1) = B32;
                        B(3, 3 * (cc)) = B33;
                        % gammaxy=du/dy+dv/dx
                        B(4, 3 * (cc) - 2) = B41;
                        B(4, 3 * (cc) - 1) = B42;
                        B(4, 3 * (cc)) = B43;
                        % gammayz=dw/dy+dv/dz
                        B(5, 3 * (cc) - 2) = B51;
                        B(5, 3 * (cc) - 1) = B52;
                        B(5, 3 * (cc)) = B53;
                        % gammaxz=dw/dx+du/dz
                        B(6, 3 * (cc) - 2) = B61;
                        B(6, 3 * (cc) - 1) = B62;
                        B(6, 3 * (cc)) = B63;
                        %%
                        Bfi(1, cc) = Nx * v11 + Ny * v12 + Nz * v13;
                        Bfi(2, cc) = Nx * v21 + Ny * v22 + Nz * v23;
                        Bfi(3, cc) = Nx * v31 + Ny * v32 + Nz * v33;
                    end
                end
            end
            cm = cm + 1;
            kufi = kufi + (B' * chargeConstantMatrixe' * Bfi) * wwwdetJ;
            kfifi = kfifi + Bfi' * permittivityMatrixEpsS * Bfi * wwwdetJ;
        end
    end
end
% ---------------------------------------------------------------------------------------------------

end
