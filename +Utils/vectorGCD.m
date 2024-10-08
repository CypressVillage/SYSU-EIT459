function myGCD = vectorGCD(inputNumbers)
% calcuate the greatest common divisor of all elements of inputNumbers
    
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
    myGCD = inputNumbers;
    return
end

if length(inputNumbers) > 2
    myGCD = Utils.vectorGCD([ inputNumbers(1:end-2), gcd( inputNumbers(end-1), inputNumbers(end)) ]);
else
    myGCD = gcd(inputNumbers(1), inputNumbers(2));
    return
end

end

