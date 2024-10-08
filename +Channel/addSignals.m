function [ sumOfSignals ] = addSignals( signal1, signal2 )
% this function adds two time domain signal that are approximately, but not
% exactly of same length. this is important when signals of two PHY
% waveforms with different roll-off factors (different filters) are added,
% for example an f-OFDM from one cell and an OFDM signal from another cell

% check inputs
if (size(signal1,1) < size(signal1,2)) || (size(signal2,1) < size(signal2,2))
    error('Input signals must be column vectors!');
end
if size(signal1,2) ~= size(signal2,2)
    if isempty(signal1)
        signal1 = zeros(size(signal2));
    elseif isempty(signal2)
        signal2 = zeros(size(signal1));
    else
        error('Number of receive antennas must match for addition!');
    end
end


% order signals by length
if size(signal1,1) == size(signal2,1)
    sumOfSignals = signal1 + signal2;
    return;
elseif size(signal1,1) > size(signal2,1)
    longSignal  = signal1;
    shortSignal = signal2;
else
    longSignal  = signal2;
    shortSignal = signal1;
end

% add zeros to short signal before adding them
nAntennas       = size(signal1,2);
diff            = size(longSignal,1) - size(shortSignal,1);
sumOfSignals    = longSignal + [zeros(floor(diff/2),nAntennas); shortSignal; zeros(ceil(diff/2),nAntennas)]; 

end

