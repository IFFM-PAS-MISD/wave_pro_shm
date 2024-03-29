function [invJ11, invJ12, invJ13, invJ21, invJ22, invJ23, invJ31, invJ32, invJ33] = ...
          inv_jacp(J11, J12, J13, J21, J22, J23, J31, J32, J33)
%
% Inverse of Jacobian in vectorized form
%
% USAGE::
%
%   [invJ11,invJ12,invJ13,invJ21,invJ22,invJ23,invJ31,invJ32,invJ33] = ...
%    inv_jacp(J11,J12,J13,J21,J22,J23,J31,J32,J33)
%
% Arguments:
%     Jij (double):
%       vectors of elements of Jacobian matrix; ij = 1,2,3
%
% Returns:
%     invJij (double):
%       inverse of Jacobian in the form of vectors; ij = 1,2,3
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------
invJ11 = (-(J23 .* J32) + J22 .* J33) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                                  (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
invJ12 = (J13 .* J32 - J12 .* J33) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                               (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
invJ13 = (-(J13 .* J22) + J12 .* J23) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                                  (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
invJ21 = (J23 .* J31 - J21 .* J33) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                               (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
invJ22 = (-(J13 .* J31) + J11 .* J33) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                                  (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
invJ23 = (J13 .* J21 - J11 .* J23) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                               (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
invJ31 = (-(J22 .* J31) + J21 .* J32) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                                  (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
invJ32 = (J12 .* J31 - J11 .* J32) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                               (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
invJ33 = (-(J12 .* J21) + J11 .* J22) ./ (-(J13 .* J22 .* J31) + (J12 .* J23 .* J31) + (J13 .* J21 .* J32) - ...
                                  (J11 .* J23 .* J32) - (J12 .* J21 .* J33) + (J11 .* J22 .* J33));
% ---------------------------------------------------------------------------------------------------

end
