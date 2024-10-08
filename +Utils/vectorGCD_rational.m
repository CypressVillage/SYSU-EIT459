function [myGCD, myGCD_numerator, myGCD_denominator] = vectorGCD_rational(numerator, denominator)
% calculate the GCD of rational numbers
% numerator and denominator must be of same length
% the element wise division numerator./denominator corresponds to the vector of input numbers

    % check input
    if ~isvector(numerator) || ~isvector(denominator)
        error('Input must be vector!');
    end

    if ~isnumeric(numerator) || ~isnumeric(denominator)
        error('Input must be numeric!');
    end
    
    if size(numerator) ~= size(denominator)
        error('Input vectors must be of same size!');
    end
    
    if any(mod(numerator,1)>0) || any(mod(denominator,1)>0)
        error('Iunput values must be integer!');
    end
    
    % find least common denominator
    comm_den = Utils.vectorLCM(denominator);
    
    % multiply all other numerators
    numerator = comm_den*numerator./denominator;
    
    % get GCD of numerator
    myGCD = Utils.vectorGCD(numerator)/comm_den;
    myGCD_numerator = Utils.vectorGCD(numerator) / gcd(Utils.vectorGCD(numerator), comm_den);
    myGCD_denominator = comm_den / gcd(Utils.vectorGCD(numerator), comm_den);
    
end

