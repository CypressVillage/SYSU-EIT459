function [ myLCM ] = vectorLCM( inputNumbers )
% this recursive funtion calculates the least common multiple of more than
% two numbers

% check input
if ~isvector(inputNumbers)
    error('Input must be vector!');
end

if ~isnumeric(inputNumbers)
    error('Input must be numeric!');
end

if any(mod(inputNumbers,1) > 0 )
    error('Input must be integers only!');
end

if isscalar(inputNumbers)
    myLCM = inputNumbers;
    return
end

if length(inputNumbers) > 2
    myLCM = Utils.vectorLCM([ inputNumbers(1:end-2), lcm( inputNumbers(end-1), inputNumbers(end)) ]);
else
    myLCM = lcm(inputNumbers(1), inputNumbers(2));
    return
end

end

