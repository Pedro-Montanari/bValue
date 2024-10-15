clear
close all
clc
%%
% Dummy parameters to create a N-element signal with amplitude between
% "a" and "b", with a final time of "tf"
N = 2000;
a = 0;
b = 4;
tf = 100; 

time= sort(0 + (tf-0).*rand(N,1));
magnitude = a + (b-a).*rand(N,1);

% plot the dummy signal
figure;
stairs(time,magnitude);
title("Original Signal");
ylabel("Signal");
xlabel("Time");

% Define the reference magnitude. The magnitude in the code is measured in
% dB, which is a relative scale. If the signal amplitude is measured in mV,
% for example, the reference amplitude is in mV. It depends on the
% characteristics of the equipment used for acquisition. 
% If the reference magnitude is not known, at least set it to be equal to
% the minimum signal amplitude, so that the magnitude will at least always 
% be positive 

% If known, for example
% args.ReferenceMagnitude = 1e-3;
% If unknown
args.ReferenceMagnitude = min(magnitude);

magnitude = log10(magnitude/args.ReferenceMagnitude);

% plot the magnitude
figure;
title("Magnitude of the signal")
stairs(time,magnitude);
ylabel("Magnitude [dB]");
xlabel("Time");

% Define b-value input arguments
args.bValueNumOfPointsInWindow = floor(N/10);

% To not utilize a sliding window set equal to 
% "args.bValueNumOfPointsInWindow"

% No sliding window
% args.bValueNumOfPointsToShiftInWindow = bValueNumOfPointsInWindow;
% Sliding window
args.bValueNumOfPointsToShiftInWindow = floor(N/10);

% To define the "bValueMagnitudeIntervalSize" it is important to take a
% look of the magnitude range of the signals. For example, if the magnitude
% levels vary between 0 and 1, a "bValueMagnitudeIntervalSize" of 0.1 could
% be used. If it varies between 0 and 2.6, for example, a 
% "bValueMagnitudeIntervalSize" of 0.5 could be used.
% To decide, one helpful way is to plot a histogram of the magnitude and to
% decide from there.

% To help deciding the "bValueMagnitudeIntervalSize"
figure;
histogram(magnitude);
title("Histogram of the signal magnitude");
ylabel("N of signals");
xlabel("Magnitude [dB]");

% After analysis of the magnitude range, for example
% args.bValueMagnitudeIntervalSize = 0.05;
% If unknown
args.bValueMagnitudeIntervalSize = 0.1;

% Run the code
b_value = bValue(time,magnitude,args);

% Plot the output data (b_value and time)
figure;
plot(b_value.Outputime,b_value.Output.Value);
title("b-Value of original signal");
ylabel("b-Value");
xlabel("Time");

% For analysis purposes, one can analyze the auxiliary data provided by the
% "Aux" substructure. It contains, for each window: 
% - The number of signals
% - The Magnitude of completeness
% - The position, in the bin, of the magnitude of completeness
%% Example 2: Signal with increasing magnitudes
% As the magnitude of signals increase, the b-value should decrease. To 
% verify if that is the case, sort the signals in ascending order
magnitude = sort(magnitude);
b_value = bValue(time,magnitude,args);

figure;
plot(b_value.Outputime,b_value.Output.Value);
title("b-Value of increasing-magnitude signal");
ylabel("b-Value");
xlabel("Time");
%% Example 3: Signal with decreasing magnitudes
% As the magnitude of signals decrease, the b-value should increase. To 
% verify if that is the case, sort the signals in descending order
magnitude = sort(magnitude,"descend");
b_value = bValue(time,magnitude,args);

figure;
plot(b_value.Outputime,b_value.Output.Value);
title("b-Value of decreasing-magnitude signal");
ylabel("b-Value");
xlabel("Time");