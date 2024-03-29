function detJ = det_jacp(J11, J12, J13, J21, J22, J23, J31, J32, J33)
%
% Determinant of Jacobian in vectorized form
%
% USAGE::
%
%   detJ = det_jacp(J11,J12,J13,J21,J22,J23,J31,J32,J33)
%
% Arguments:
%     Jij (double):
%       vectors of elements of Jacobian matrix; ij = 1,2,3
%
% Returns:
%     detJ (double):
%       determinant of Jacobian in the form of vector
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------
detJ = (J11 .* J22 .* J33 + J12 .* J23 .* J31 + J13 .* J21 .* J32 - J13 .* J22 .* J31 - J11 .* ...
        J23 .* J32 - J12 .* J21 .* J33);
% ---------------------------------------------------------------------------------------------------

end
