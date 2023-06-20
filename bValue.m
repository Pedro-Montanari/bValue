function b_value = bValue(T,args)
%bValue : Calculates the b-value and associated parameters
%   Updated in 06/20/2023
%   This function calculates the b-value for a set of data.
%   Input arguments:
%       T : (structure) containing the following field:
%           1.  Time (double) = array containing signal times
%       args = structure containing the following fields:
%           1. MagRange (double) : scalar containing the difference between
%           maximum and minimum
%           signal magnitude
%           2. Magnitude (double) : array containing signal magnitude
%           3. bValueNumOfPointsInWindow (double) : scalar indicating how
%           many values to use in each b-value calculation
%           4. bValueNumOfPointsToShiftInWindow (double) : scalar
%           indicating how many points to move each window.
%           The calculation uses a sliding window scheme, which means
%           overlap can be taken into account. For no overlap specify
%           "bValueNumOfPointsToShiftInWindow" =
%           "bValueNumOfPointsInWindow"
%   Output arguments:
%       b_value : (structure) with the following fields:
%           1. Value : b_value itself
%           2. Magnitude : Magnitude of signals (optional)
%           3. NofHitsInEachInterval : Number of AE hits in each magnitude
%           interval defined (optional)
%           4. MagnitudeOfCompleteness : Magnitude of Completeness of the
%           signals (see 10.1016/j.chaos.2015.09.004 for info) (optional)
%           5. PositionOfMagnitudeOfCompleteness : Position of the
%           magnitude which is the magnitude of completeness (optional)

b_value(1).Time(:,1) = b_value_time_points;

b_value(1).Magnitude = CurrentWindowMagnitude;

b_value(1).NofHitsInEachInterval = NofHitsInEachInterval;

b_value(1).MagnitudeOfCompleteness = MagnitudeOfCompleteness;

b_value(1).PositionOfMagnitudeOfCompleteness = PositionOfMagnitudeOfCompleteness;
%   NOTE: "T" is a structure because originally this function was designed
%   to work with a code where it was a structure, but for this function
%   actually only a time array "T" is necessary.

b_value = struct('Value',{},'Time',{},'Magnitude',{},'NofHitsInEachInterval',{},'MagnitudeOfCompleteness',{},'PositionOfMagnitudeOfCompleteness',{});

NumOfIntervals = floor(args.MagRange/args.bValueMagnitudeIntervalSize);

bValueIntervals = linspace(min(args.Magnitude),max(args.Magnitude),NumOfIntervals);

NofHitsInEachInterval = cell(1,NumOfIntervals);

CurrentWindowPosition = 1:args.bValueNumOfPointsInWindow;

counter = 1;

while(CurrentWindowPosition(end) <= numel(T.Time))
    NofHitsInEachInterval{counter} = zeros(2,NumOfIntervals);

    for iIntervals=1:NumOfIntervals
        NofHitsInEachInterval{counter}(1,iIntervals) = bValueIntervals(iIntervals);

        NofHitsInEachInterval{counter}(2,iIntervals) = nnz(args.Magnitude(CurrentWindowPosition) >= bValueIntervals(iIntervals));
    end

    CurrentWindowMagnitude{counter} = args.Magnitude(CurrentWindowPosition);

    b_value_time_points(counter) = mean(T.Time(CurrentWindowPosition));

    [NofHitsOfMagnitudeOfCompleteness(counter),PositionOfMagnitudeOfCompleteness(counter)] = max(flip(NofHitsInEachInterval{counter}(2,:)));

    PositionOfMagnitudeOfCompleteness(counter) = numel(NofHitsInEachInterval{counter}(2,:)) - PositionOfMagnitudeOfCompleteness(counter) + 1;

    MagnitudeOfCompleteness(counter) = NofHitsInEachInterval{counter}(1,PositionOfMagnitudeOfCompleteness(counter));

    CurrentWindowPosition = CurrentWindowPosition + args.bValueNumOfPointsToShiftInWindow;

    counter = counter+1;
end

RegressionLineCoefficients = zeros(numel(b_value_time_points),2);

X = cell(1,numel(b_value_time_points));

Y = cell(1,numel(b_value_time_points));

yl = cell(1,numel(b_value_time_points));

for iBvalue=1:numel(b_value_time_points)
    Filt = NofHitsInEachInterval{iBvalue}(2,:) ~= 0;

    X{iBvalue} = NofHitsInEachInterval{iBvalue}(1,Filt);

    Y{iBvalue} = log10(NofHitsInEachInterval{iBvalue}(2,Filt));

    RegressionLineCoefficients(iBvalue,:) = polyfix(X{iBvalue},Y{iBvalue},1,X{iBvalue}(PositionOfMagnitudeOfCompleteness(iBvalue)), ...
        ...
        Y{iBvalue}(PositionOfMagnitudeOfCompleteness(iBvalue)));

    xl{iBvalue} = X{iBvalue}(PositionOfMagnitudeOfCompleteness(iBvalue):end);

    yl{iBvalue} = RegressionLineCoefficients(iBvalue,1)*xl{iBvalue} + RegressionLineCoefficients(iBvalue,2);
end

b_value(1).Value(:,1) = abs(RegressionLineCoefficients(:,1));

b_value(1).Time(:,1) = b_value_time_points;

b_value(1).Magnitude = CurrentWindowMagnitude;

b_value(1).NofHitsInEachInterval = NofHitsInEachInterval;

b_value(1).MagnitudeOfCompleteness = MagnitudeOfCompleteness;

b_value(1).PositionOfMagnitudeOfCompleteness = PositionOfMagnitudeOfCompleteness;
end