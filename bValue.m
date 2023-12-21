function b_value = bValue(T,args)
%bValue : Calculates the b-value and associated parameters
%   This function calculates the b-value for a set of data.
%   Input arguments:
%       T : (structure) containing the following field:
%           1.  Time (double) = array containing signal times
%           2.  Magnitude(double) = array containing the signal Magnitude
%           in dB
%       args = structure containing the following fields:
%           1. bValueNumOfPointsInWindow (double) : scalar indicating how
%           many values to use in each b-value calculation
%           2. bValueNumOfPointsToShiftInWindow (double) : scalar
%           indicating how many points to move each window.
%           The calculation uses a sliding window scheme, which means
%           overlap can be taken into account. For no overlap specify
%           'bValueNumOfPointsToShiftInWindow' =
%           'bValueNumOfPointsInWindow'
%           3.bValueMagnitudeIntervalSize (double): scalar indicating the
%           interval size of magnitudes. By default use 0.1.
%   Output arguments:
%       b_value : (structure) with the following fields:
%           1. Input : (structure) with the following fields
%               i.  Time: (double) array of the time values of the 
%                   original signal
%               ii. Signal: (double) array of the original input signal
%               iii. Magnitude: (double) array of the magnitudes of the
%               original signal
%           2. Output: (structure) with the following fields
%               i. Value : b_value itself
%               ii.Time: time array of the b_values (obs.: not a "real"
%               time, it is the time in the middle of each window.
%           3. Aux: (structure) with the following fiels
%               i. NofHitsInEachInterval : Number of AE hits in each magnitude
%           interval defined
%               ii. MagnitudeOfCompleteness : Magnitude of Completeness of the
%           signals (see 10.1016/j.chaos.2015.09.004 for info)
%               iii. PositionOfMagnitudeOfCompleteness : Position of the
%           magnitude which is the magnitude of completeness

args.MagRange = range(T.Magnitude);

NumOfIntervals = floor(args.MagRange/args.bValueMagnitudeIntervalSize);
bValueIntervals(:,1) = ...
    linspace(min(T.Magnitude),max(T.Magnitude),NumOfIntervals);
NofHitsInEachInterval = cell(NumOfIntervals,1);
CurrentWindowPosition = 1:args.bValueNumOfPointsInWindow;
counter = 1;

% Initialize variables
CurrentWindowMagnitudes = cell(numel(T.Signal),1);
bValueTimePoints = zeros(numel(T.Signal),1);
NofHitsOfMagnitudeOfCompleteness = zeros(numel(T.Signal),1);
PositionOfMagnitudeOfCompleteness = zeros(numel(T.Signal),1);
MagnitudeOfCompleteness = zeros(numel(T.Signal),1);

while(CurrentWindowPosition(end) <= numel(T.Time))
    NofHitsInEachInterval{counter} = zeros(2,NumOfIntervals);

    for iIntervals=1:NumOfIntervals
        NofHitsInEachInterval{counter}(1,iIntervals) = ...
            bValueIntervals(iIntervals);
        NofHitsInEachInterval{counter}(2,iIntervals) = ...
            nnz(T.Magnitude(CurrentWindowPosition) >= bValueIntervals(iIntervals));
    end

    CurrentWindowMagnitudes{counter} = T.Magnitude(CurrentWindowPosition);
    bValueTimePoints(counter) = mean(T.Time(CurrentWindowPosition));
    
    [NofHitsOfMagnitudeOfCompleteness(counter),...
     PositionOfMagnitudeOfCompleteness(counter)] = ...
        max(flip(NofHitsInEachInterval{counter}(2,:)));

    PositionOfMagnitudeOfCompleteness(counter) = ...
        numel(NofHitsInEachInterval{counter}(2,:)) - ...
        PositionOfMagnitudeOfCompleteness(counter) + 1;
    MagnitudeOfCompleteness(counter) = ...
        NofHitsInEachInterval{counter}(1,PositionOfMagnitudeOfCompleteness(counter));
    CurrentWindowPosition = ...
        CurrentWindowPosition + args.bValueNumOfPointsToShiftInWindow;
    counter = counter+1;
end

% Remove empty data
NofHitsInEachInterval = ...
    NofHitsInEachInterval(~cellfun(@isempty,NofHitsInEachInterval));
bValueTimePoints(counter:end) = [];
MagnitudeOfCompleteness(counter:end) = [];
PositionOfMagnitudeOfCompleteness(counter:end) = [];

RegressionLineCoefficients = zeros(numel(bValueTimePoints),2);
X = cell(numel(bValueTimePoints),1);
Y = cell(numel(bValueTimePoints),1);
xl = cell(numel(bValueTimePoints),1);
yl = cell(numel(bValueTimePoints),1);

for iBvalue=1:numel(bValueTimePoints)
    Filt = NofHitsInEachInterval{iBvalue}(2,:) ~= 0;
    X{iBvalue} = NofHitsInEachInterval{iBvalue}(1,Filt);
    Y{iBvalue} = log10(NofHitsInEachInterval{iBvalue}(2,Filt));
    RegressionLineCoefficients(iBvalue,:) = ...
        polyfix(X{iBvalue},...
                Y{iBvalue},1,...
                X{iBvalue}(PositionOfMagnitudeOfCompleteness(iBvalue)), ...
                Y{iBvalue}(PositionOfMagnitudeOfCompleteness(iBvalue)));
    xl{iBvalue} = X{iBvalue}(PositionOfMagnitudeOfCompleteness(iBvalue):end);
    yl{iBvalue} = RegressionLineCoefficients(iBvalue,1)*xl{iBvalue} + ...
                  RegressionLineCoefficients(iBvalue,2);
end

b_value.Input.Time = T.Time;
b_value.Input.Signal = T.Signal;
b_value.Input.Magnitude = T.Magnitude;
b_value.Output.Value(:,1) = abs(RegressionLineCoefficients(:,1));
b_value.Output.Time(:,1) = bValueTimePoints;
b_value.Aux.NofHitsInEachInterval = cellfun(@transpose,...
                                           NofHitsInEachInterval,...
                                           'UniformOutput', false);
b_value.Aux.MagnitudeOfCompleteness = MagnitudeOfCompleteness;
b_value.Aux.PositionOfMagnitudeOfCompleteness = ...
    PositionOfMagnitudeOfCompleteness;
end