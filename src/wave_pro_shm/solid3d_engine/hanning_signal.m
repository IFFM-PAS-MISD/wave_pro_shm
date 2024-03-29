function [t, st] = hanning_signal(dt, nSamples, modulationFrequency, carrierFrequency, initiationTime)
%
% Calculate Hann windowed signal
%
% USAGE::
%
%   [t,st] = hanning_signal(dt,nSamples,modulationFrequency,carrierFrequency,initiationTime)
%
% Arguments:
%     dt (double):
%       time step, Units [s]
%
%     nSamples (integer):
%       number of samples
%
%     modulationFrequency (double):
%       modulation frequency, Units [Hz]
%
%     carrierFrequency (double):
%       carrier frequency, Units [Hz]
%
%     initiationTime (double):
%       initiation time (delay), Units [s]
%
% Returns:
%     t (double):
%       time vector, Units [s], dimensions [nSamples,1]
%
%     st (double):
%       time domain excitation signal in range [-1,1], dimensions [nSamples,1]
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

t_t = 1 / modulationFrequency;   % total duration time of the excitation [s]
t_2 = initiationTime + t_t;      % excitation termination time [s]

st = zeros(nSamples, 1);
t = zeros(nSamples, 1);
for i = 1:nSamples
    st(i) = 0.0;
    t(i) = (i - 1) * dt;
    if (t(i) >= initiationTime) && (t(i) <= t_2)
        st(i) = 0.5 * (1 - cos(2 * pi * modulationFrequency * (t(i) - t_2))) * ...
                       sin(2 * pi * carrierFrequency * (t(i) - initiationTime));
    end
end

% ---------------------------------------------------------------------------------------------------

end
