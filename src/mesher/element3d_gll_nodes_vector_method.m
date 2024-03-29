function [Xe, Ye, Ze] = element3d_gll_nodes_vector_method(ksi, eta, dzeta, x, y, z)
%
% Coordinates of the nodes of 3D spectral element according to Gauss-Lobatto-Legendre nodes
% nodes distribution is specified by Gauss-Lobatto-Legendre (GLL) points
% and can be different along each axis (ksi, eta, dzeta)
% procedure valid only for linear geometry
%
% USAGE::
%
%   [Xe, Ye, Ze] = element3d_gll_nodes_vector_method(ksi, eta, dzeta, x, y, z)
%
% Arguments:
%     ksi (double):
%       GLL distribution along local element x axis
%       range -1:1, dimensions [1, nElementNodesX]
%
%     eta (double):
%       GLL distribution along local element y axis
%       range -1:1, dimensions [1, nElementNodesY]
%
%     dzeta (double):
%       GLL distribution along local element z axis
%       range -1:1, dimensions [1, nElementNodesZ]
%
%     x,y,z (double):
%      coordinates of corner nodes of the element, dimensions [1, 8], Units [m]
%
% Returns:
%     Xe (double):
%       global x-axis coordinates of all nodes of the element,
%       dimensions [1, nElementNodesX*nElementNodesY*nElementNodesZ], Units [m]
%
%     Ye (double):
%       global y-axis coordinates of all nodes of the element,
%       dimensions [1, nElementNodesX*nElementNodesY*nElementNodesZ], Units [m]
%
%     Ze (double):
%       global z-axis coordinates of all nodes of the element,
%       dimensions [1, nElementNodesX*nElementNodesY*nElementNodesZ], Units [m]
%
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

nElementNodesX = length(ksi);
nElementNodesY = length(eta);
nElementNodesZ = length(dzeta);

Xe = zeros(1, nElementNodesX * nElementNodesY * nElementNodesZ);
Ye = zeros(1, nElementNodesX * nElementNodesY * nElementNodesZ);
Ze = zeros(1, nElementNodesX * nElementNodesY * nElementNodesZ);

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
x5 = x(5);
x6 = x(6);
x7 = x(7);
x8 = x(8);
y1 = y(1);
y2 = y(2);
y3 = y(3);
y4 = y(4);
y5 = y(5);
y6 = y(6);
y7 = y(7);
y8 = y(8);
z1 = z(1);
z2 = z(2);
z3 = z(3);
z4 = z(4);
z5 = z(5);
z6 = z(6);
z7 = z(7);
z8 = z(8);

% vectors along element z axis
% edge 1-5
h1x = x5 - x1;
h1y = y5 - y1;
h1z = z5 - z1;
% edge 2-6
h2x = x6 - x2;
h2y = y6 - y2;
h2z = z6 - z2;
% edge 3-7
h3x = x7 - x3;
h3y = y7 - y3;
h3z = z7 - z3;
% edge 4-8
h4x = x8 - x4;
h4y = y8 - y4;
h4z = z8 - z4;

counter = 0;
for iNodeZ = 1:nElementNodesZ % z axis
    % gll coordinates at edge 1-5
    X1_ = x1 + h1x * (1 + dzeta(iNodeZ)) / 2;
    Y1_ = y1 + h1y * (1 + dzeta(iNodeZ)) / 2;
    Z1_ = z1 + h1z * (1 + dzeta(iNodeZ)) / 2;

    % coordinates at edge 2-6
    X2_ = x2 + h2x * (1 + dzeta(iNodeZ)) / 2;
    Y2_ = y2 + h2y * (1 + dzeta(iNodeZ)) / 2;
    Z2_ = z2 + h2z * (1 + dzeta(iNodeZ)) / 2;

    % coordinates at edge 3-7
    X3_ = x3 + h3x * (1 + dzeta(iNodeZ)) / 2;
    Y3_ = y3 + h3y * (1 + dzeta(iNodeZ)) / 2;
    Z3_ = z3 + h3z * (1 + dzeta(iNodeZ)) / 2;

    % coordinates at edge 4-8
    X4_ = x4 + h4x * (1 + dzeta(iNodeZ)) / 2;
    Y4_ = y4 + h4y * (1 + dzeta(iNodeZ)) / 2;
    Z4_ = z4 + h4z * (1 + dzeta(iNodeZ)) / 2;

    % vectors along element y axis
    % edge 1-3
    bx = X3_ - X1_;
    by = Y3_ - Y1_;
    bz = Z3_ - Z1_;
    % edge 2-4
    cx = X4_ - X2_;
    cy = Y4_ - Y2_;
    cz = Z4_ - Z2_;
    for jNodeY = 1:nElementNodesY % y axis
        % gll coordinates along y axis
        Xb_ = X1_ + bx * (1 + eta(jNodeY)) / 2;
        Yb_ = Y1_ + by * (1 + eta(jNodeY)) / 2;
        Zb_ = Z1_ + bz * (1 + eta(jNodeY)) / 2;

        Xc_ = X2_ + cx * (1 + eta(jNodeY)) / 2;
        Yc_ = Y2_ + cy * (1 + eta(jNodeY)) / 2;
        Zc_ = Z2_ + cz * (1 + eta(jNodeY)) / 2;

        % vectors along x axis
        ax = Xc_ - Xb_;
        ay = Yc_ - Yb_;
        az = Zc_ - Zb_;
        for kNodeX = 1:nElementNodesX % x axis
            counter = counter + 1;
            % gll coordinates along x axis
            Xe(counter) = Xb_ + ax * (1 + ksi(kNodeX)) / 2;
            Ye(counter) = Yb_ + ay * (1 + ksi(kNodeX)) / 2;
            Ze(counter) = Zb_ + az * (1 + ksi(kNodeX)) / 2;
        end
    end
end
% ---------------------------------------------------------------------------------------------------

end
