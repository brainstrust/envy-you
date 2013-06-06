function [response] = createSig(Time, SwitchingT, LowerBound, UpperBound, Ramp)
% This Mfile creates a smooth ramp function output for a given time step or 
% time vector based on the sigmodal behavior of a hyperbolic tangens.
%
% Default values SwitchingX = 0, LowerBound = 0, UpperBound=1, Ramp=3
%
% LowerBound < UpperBound Behavior = ,,,,/'''''
% UpperBound < LowerBound Behavior = ''''\,,,,,
%
% Time response will be slower by incresing the Ramp values, and wil be
% quicker by decreasing the values towards zero. Hence that an overshoot
% can be created when values are not matched nicely.
%
% Evert de Kock 2012

%% Check input arguments and set to defaults if necessary
if nargin <= 1; SwitchingT = 0;end
if nargin <= 2; LowerBound = 0;end
if nargin <= 3; UpperBound = 1;end
if nargin <= 4; Ramp       = 3;end

%% VectorInputCheck
if isscalar(Time)
    vectorSize = 1;
elseif isvector(Time)
        vectorSize = ones(length(Time),1);
else
    fprintf('Warning: Input parameter Time (size [%i, %i]) in createRef exceedes vector dimensions\n',size(Time,1), size(Time,2))
    % take the first element of time to continue calculations
    Time = Time(1);
    vectorSize = 1;
end

%% Actual sigmodial behavior
response  = ((UpperBound - LowerBound)/(2))...               % Half amplitude
             *(1+tanh((Time-SwitchingT*vectorSize)/Ramp))... % actual bimodal behavior
             +LowerBound*vectorSize;                         % ofsett wrt the x axis
end

