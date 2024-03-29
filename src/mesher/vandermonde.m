function [vandermondeMatrix] = vandermonde(ksi)
%
% Calculate Vandermonde matrix
%
% USAGE::
%
%   [vandermondeMatrix] = vandermonde(ksi)
%
% Arguments:
%     ksi (double):
%       local coordinates in range -1:1, row vector of dimensions [1,n]
%       where n is the number of nodes in 1D spectral element
%
% Returns:
%     vandermondeMatrix (double):
%       Vandermonde matrix, dimensions [n,n]
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

vandermondeMatrix = zeros(length(ksi), length(ksi));

for iRow = 1:length(ksi)
    for jCol = 1:length(ksi)
        vandermondeMatrix(iRow, jCol) = ksi(jCol)^(iRow - 1);
    end
end

% ---------------------------------------------------------------------------------------------------

end
