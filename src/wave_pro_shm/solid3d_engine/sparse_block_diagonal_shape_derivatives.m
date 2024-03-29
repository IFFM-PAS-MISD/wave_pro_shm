function [Nksi_P, Neta_P, Ndzeta_P] = sparse_block_diagonal_shape_derivatives(shapePrimKsi_e, ...
                                                                              shapePrimEta_e, ...
                                                                              shapePrimDzeta_e, ...
                                                                              nElements, ...
                                                                              nElement3DNodes)
%
% create sparse block-diagonal matrix of shape function derivatives
%
% USAGE::
%
%   [Nksi_P,Neta_P,Ndzeta_P] = sparse_block_diagonal_shape_derivatives(shapePrimKsi_e,...
%                              shapePrimEta_e,shapePrimDzeta_e,nElements,nElement3DNodes)
%
% Arguments:
%     shapePrimKsi_e (double):
%       derivatives of shape functions in respect to ksi (nodal values)
%
%     shapePrimEta_e (double):
%       derivatives of shape functions in respect to eta (nodal values)
%
%     shapePrimDzeta_e (double):
%       derivatives of shape functions in respect to dzeta
%
%     nElements (integer):
%       number of elements
%
%     nElement3DNodes (integer):
%       number of nodes in a 3D solid spectral element
%
% .. Note::
%    _e at the end means that this is for a single element
%
% Returns:
%     Nksi_P (double):
%       sparse block-diagonal matrix of shape function derivatives along ksi axis
%
%     Neta_P (double):
%       sparse block-diagonal matrix of shape function derivatives along ksi axis
%
%     Ndzeta_P (double):
%       sparse block-diagonal matrix of shape function derivatives along dzeta axis
%
%
% .. Note::
%    _P at the end means that it in vectorized form intended for parallel computation
%
%
% .. seealso:: P. Kudela, Parallel implementation of spectral element method for Lamb wave
%              propagation modeling, *International Journal for Numerical Methods in Engineering*,
%              2016, 106:413-429, `doi:10.1002/nme.5119 <http://doi.wiley.com/10.1002/nme.5119>`_ ,
%              Eqs. (2.2)-(2.4)
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

[indxi, indxj] = find(shapePrimKsi_e); % indices of nonzero values for derivative in ksi direction
[indyi, indyj] = find(shapePrimEta_e); % indices of nonzero values for derivative in eta direction
[indzi, indzj] = find(shapePrimDzeta_e); % indices of nonzero values for derivative in dzeta direction

ixl = length(indxi);
iyl = length(indyi);
izl = length(indzi);
Nksi_Ps = zeros(ixl, 1);
Neta_Ps = zeros(iyl, 1);
Ndzeta_Ps = zeros(izl, 1);
for k = 1:ixl
    Nksi_Ps(k) = shapePrimKsi_e(indxi(k), indxj(k));
end
for k = 1:iyl
    Neta_Ps(k) = shapePrimEta_e(indyi(k), indyj(k));
end
for k = 1:izl
    Ndzeta_Ps(k) = shapePrimDzeta_e(indzi(k), indzj(k));
end
Nxs = repmat(Nksi_Ps, nElements, 1);
Nys = repmat(Neta_Ps, nElements, 1);
Nzs = repmat(Ndzeta_Ps, nElements, 1);
clear Nksi_Ps Neta_Ps Ndzeta_Ps;
Indxi = zeros(nElements * ixl, 1);
Indxj = zeros(nElements * ixl, 1);
Indyi = zeros(nElements * iyl, 1);
Indyj = zeros(nElements * iyl, 1);
Indzi = zeros(nElements * izl, 1);
Indzj = zeros(nElements * izl, 1);
for iElement = 1:nElements
    off = (iElement - 1) * nElement3DNodes;
    n1 = (iElement - 1) * ixl + 1;
    n2 = n1 + ixl - 1;
    n3 = (iElement - 1) * iyl + 1;
    n4 = n3 + iyl - 1;
    n5 = (iElement - 1) * izl + 1;
    n6 = n5 + izl - 1;
    Indxi(n1:n2, 1) = off + indxi;
    Indxj(n1:n2, 1) = off + indxj;
    Indyi(n3:n4, 1) = off + indyi;
    Indyj(n3:n4, 1) = off + indyj;
    Indzi(n5:n6, 1) = off + indzi;
    Indzj(n5:n6, 1) = off + indzj;
end
% abbreviate naming: Nksi is shape function derivative in respect to ksi
Nksi_P = sparse(Indxi, Indxj, Nxs);
clear Nxs;
Neta_P = sparse(Indyi, Indyj, Nys);
clear Nys;
Ndzeta_P = sparse(Indzi, Indzj, Nzs);
clear Nzs;
clear indxi indxj indyi indyj indzi indzj;
% ---------------------------------------------------------------------------------------------------
end
