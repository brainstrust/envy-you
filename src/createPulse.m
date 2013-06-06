function [response] = createPulse(Time, SwitchingLeft, SwitchingRight, LowerBound, UpperBound, RampLeft, RampRight)
% This Mfile creates a smooth pulse function output for a given time step or 
% time vector based on the sigmodal behavior of two hyperbolic tangens.
%
% Default values: SwitchingLeft = 0
%                 SwitchingRight = 1
%                 LowerBound    = 0
%                 UpperBound    = 1
%                 RampLeft      =(1/3)
%                 RampRight     =(1/3)
%
% Time response will be slower by incresing the Ramp values, and wil be
% quicker by decreasing the values towards zero. Hence that an overshoot
% can be created when values are not matched nicely.
%
% LowerBound < UpperBound Behavior = ,,,,/'''''\,,,,
% UpperBound < LowerBound Behavior = ''''\,,,,,/''''
%
% Evert de Kock 2012
%% Check input arguments and set to defaults if necessary
if nargin <= 1; SwitchingLeft   = 0;    end
if nargin <= 2; SwitchingRight  = 1;    end
if nargin <= 3; LowerBound      = 0;    end
if nargin <= 4; UpperBound      = 1;    end
if nargin <= 5; RampLeft        = 3;    end
if nargin <= 6; RampRight       = 3;    end

%% VectorInputCheck
if isscalar(Time)
    vectorSize = 1;
elseif isvector(Time)
    vectorSize = ones(1,length(Time));
else
    fprintf('Warning: Input parameter Time (size [%i, %i]) in createRef exceedes vector dimensions\n',size(Time,1), size(Time,2))
    % take the first element of time to continue calculations
    Time = Time(1);
    vectorSize = 1;
end

%% Check SwitchLeft < SwitchRight
if SwitchingLeft > SwitchingRight
    dumpVar = SwitchingLeft;
    SwitchingLeft = SwitchingRight;
    SwitchingRight = dumpVar;
    fprintf('Warning: Left must be smaller than Right, Copy That? No Worries, it is fixxed it for you ;-) \n')
end
%% Actual Pulse Calculation
% Calculate first increasing part of the pulse:
ResLeft  = createSig(Time', SwitchingLeft,  LowerBound, UpperBound, RampLeft);
% Calculate second decreasing part of the pulse:
ResRight = createSig(Time', SwitchingRight, UpperBound, LowerBound, RampRight);

% plot(Time,ResLeft,'b*')
% plot(Time,ResRight,'g*')

% Add both suppulses to one pulse: [[Little bit dirty programming,,, might be fixxed]]
if LowerBound < UpperBound
response = (ResLeft + ResRight) - vectorSize'*UpperBound ;
elseif LowerBound > UpperBound
response = (ResLeft + ResRight) - vectorSize'*UpperBound ;
else
fprintf('Warning: Pulse function boundries are equal no effective pules created\n')
response = vectorSize'*Lowerbound;
end
end