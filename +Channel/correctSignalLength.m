function [ correctedSignal ] = correctSignalLength( signal, signalLength )
% this function shapes the signal to have aspecified length in a symmetrix
% way

% check input
if (size(signal,1) < size(signal,2))
    error('Input signal must be column vector!');
end

nAntennas = size(signal,2);
% correct signal length
if size(signal,1) == signalLength
    correctedSignal = signal;
elseif size(signal,1) > signalLength
    diff = size(signal,1) - signalLength;
    correctedSignal = signal(floor(diff/2)+1:floor(diff/2)+signalLength,:);
else
    diff = signalLength - size(signal,1);
    correctedSignal = [zeros(floor(diff/2),nAntennas); signal; zeros(ceil(diff/2),nAntennas)];
end


end

