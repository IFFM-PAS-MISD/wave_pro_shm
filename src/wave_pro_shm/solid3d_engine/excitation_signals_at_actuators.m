function Excitation = excitation_signals_at_actuators(Excitation)
%
% Calculate excitation signals corresponding to each actuator
%
% USAGE::
%
%   Excitation = excitation_signals(Excitation)
%
% Arguments:
%     Excitation (struct):
%       structure describing excitation with fields
%
%       assignment (obj):
%         excitation assignment to nActuators with fields
%
%         category (integer):
%           1-Hann, 2-Chirp, 3-Gauss
%
%         index (integer):
%           signal index corresponding to iActuator
%
%       nSamples (integer):
%         number of samples
%
%       nActuators (integer):
%         number of actuators
%
%       totalTime (double):
%         total wave propagation time, Units [s]
%
%       Hann (obj):
%         Hann signal objects with category specific fields
%
%         peakVoltage (double):
%           signal peak voltage, Units [V]
%
%         modulationFrequency (double):
%           vector of modulation frequencies for signals at nActuators, Units [Hz]
%
%         carrierFrequency (double):
%           vector of carrier frequencies for signals at nActuators Units [s]
%
%         initiationTime (double):
%           vector of initiation times (delays) for signals at nActuators Units [Hz]
%
% Returns:
%     Excitation (struct):
%       structure expanded by fields
%
%         signals (double): time domain excitation signals in range [-1,1],
%                           dimensions [nSamples,nActuators]
%
%         timeVector (double):
%           time vector, Units [s]
%
%         dt (double):
%           time step, Units [s]
%
% .. seealso:: Function :func:`hanning_signal`
%
% TODO: provisions for Chirp and Gauss signals are made which should be implemented
%
%
% (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl
% Institute of Fluid Flow Machinery, Polish Academy of Sciences
% Mechanics of Intelligent Structures Department

% ---------------------------------------------------------------------------------------------------

Excitation.signals = zeros(Excitation.nSamples, Excitation.nActuators);

for iActuator = 1:Excitation.nActuators
    if Excitation.assignment(iActuator).category == 1  % Hann signal
        iHannSignal = Excitation.assignment(iActuator).index;
        if Excitation.Hann(iHannSignal).peakVoltage == 0
            signal = zeros(Excitation.nSamples, 1);
        else
            [~, signal] = hanning_signal(Excitation.dt, Excitation.nSamples, ...
                                        Excitation.Hann(iHannSignal).modulationFrequency, ...
                                        Excitation.Hann(iHannSignal).carrierFrequency, ...
                                        Excitation.Hann(iHannSignal).initiationTime);
        end

    elseif Excitation.assignment(iActuator).category == 2  % Chirp signal
        iChirpSignal = Excitation.assignment(iActuator).index;
    elseif Excitation.assignment(iActuator).category == 3  % Gauss signal
        iGaussSignal = Excitation.assignment(iActuator).index;
    end
    Excitation.signals(:, iActuator) = signal(:, 1);
end

timeVector = 0:Excitation.dt:(Excitation.nSamples - 1) * Excitation.dt;
Excitation.timeVector = timeVector;
% ---------------------------------------------------------------------------------------------------

end
