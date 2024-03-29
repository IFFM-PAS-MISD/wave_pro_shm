function [xNew, yNew, zNew] = nodal_base_change_3d(ksiNew, etaNew, dzetaNew, ...
                                                 ksiOld, etaOld, dzetaOld, ...
                                                 xOld, yOld, zOld)
%
% Change nodal basis of 3D solid element (node coordinates)
% Applicable for elements geometry approximated by any high-order Lagrange function
%
% USAGE::
%
%   [xNew,yNew,zNew] = nodal_base_change_3D(ksiNew,etaNew,dzetaNew,...
%                                           ksiOld,etaOld,dzetaOld,...
%                                           xOld,yOld,zOld)
%
% Arguments:
%     ksiNew (double):
%       new nodal basis along ksi, coordinates in range -1:1
%
%     etaNew (double):
%       new nodal basis along eta, coordinates in range -1:1
%
%     dzetaNew (double):
%       new nodal basis along dzeta, coordinates in range -1:1
%
%     ksiOld (double):
%       old nodal basis along ksi, coordinates in range -1:1
%
%     etaOld (double):
%       old nodal basis along eta, coordinates in range -1:1
%
%     dzetaOld (double):
%       old nodal basis along dzeta, coordinates in range -1:1
%
%     xOld, yOld, zOld (double):
%       global coordinates in old nodal basis
%
% Returns:
%     xNew, yNew, zNew (double):
%       global coordinates in new nodal basis
%
% TODO: vectorize code
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

[vandermondeMatrixOldX] = vandermonde(ksiOld);
[vandermondeMatrixOldY] = vandermonde(etaOld);
[vandermondeMatrixOldZ] = vandermonde(dzetaOld);

inverseVandermondeMatrixOldX = inv(vandermondeMatrixOldX);
inverseVandermondeMatrixOldY = inv(vandermondeMatrixOldY);
inverseVandermondeMatrixOldZ = inv(vandermondeMatrixOldZ);

[vandermondeMatrixNewX] = vandermonde(ksiNew);
[vandermondeMatrixNewY] = vandermonde(etaNew);
[vandermondeMatrixNewZ] = vandermonde(dzetaNew);

NodalBaseNewX = inverseVandermondeMatrixOldX * vandermondeMatrixNewX;
NodalBaseNewY = inverseVandermondeMatrixOldY * vandermondeMatrixNewY;
NodalBaseNewZ = inverseVandermondeMatrixOldZ * vandermondeMatrixNewZ;

xNew = zeros(length(ksiNew) * length(etaNew) * length(dzetaNew), 1);
yNew = zeros(length(ksiNew) * length(etaNew) * length(dzetaNew), 1);
zNew = zeros(length(ksiNew) * length(etaNew) * length(dzetaNew), 1);

cc = 0;
for i1 = 1:length(dzetaOld)
    for j1 = 1:length(etaOld)
        for k1 = 1:length(ksiOld)
            cc = cc + 1;
            c = 0;
            for i2 = 1:length(dzetaNew)
                for j2 = 1:length(etaNew)
                    for k2 = 1:length(ksiNew)
                        c = c + 1;
                        xNew(c) = xNew(c) + ...
                            NodalBaseNewZ(i1, i2) * NodalBaseNewY(j1, j2) * NodalBaseNewX(k1, k2) * xOld(cc);
                        yNew(c) = yNew(c) + ...
                            NodalBaseNewZ(i1, i2) * NodalBaseNewY(j1, j2) * NodalBaseNewX(k1, k2) * yOld(cc);
                        zNew(c) = zNew(c) + ...
                            NodalBaseNewZ(i1, i2) * NodalBaseNewY(j1, j2) * NodalBaseNewX(k1, k2) * zOld(cc);
                    end
                end
            end
        end
    end
end

% ---------------------------------------------------------------------------------------------------

end
