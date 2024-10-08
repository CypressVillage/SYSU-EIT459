classdef CRC < handle
% Cyclic Redundancy Check Class (CRC)
%
%   Author
%       - Bashar Tahir, btahir@nt.tuwien.ac.at
%         (c) 2018 Institute of Telecommunications, TU Wien.
%         www.nt.tuwien.ac.at 
% 
%   Performs CRC calclation and attachment.
%   Performs CRC error detection.
% 
%   The class supports CRC polynomials based on LTE/5G-NR.
 
properties (SetAccess = private, Hidden = true)
    % CRC ploymoials from LTE/5G-NR
    
    CRC24A = [1,1,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,1,1,1,1,1,0,1,1];
    CRC24B = [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1];
    CRC24C = [1,1,0,1,1,0,0,1,0,1,0,1,1,0,0,0,1,0,0,0,1,0,1,1,1];
    CRC16 = [1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1];
    CRC11 = [1,1,1,0,0,0,1,0,0,0,0,1];
    CRC8 = [1,1,0,0,1,1,0,1,1];
end

methods (Static)
%%
    function CRCAttached = Attach(varargin)
    % Calculates the CRC and attaches it to the input block
    % Calculates the CRC and attaches it to the input block 
    % based on the selected polynomial.
        obj = Coding.ChannelCodingComponents.CRC;
        InputBits = varargin{1};
        switch varargin{2}
            case '24A'
                CRCPoly = obj.CRC24A;
            case '24B'
                CRCPoly = obj.CRC24B;
            case '24C'
                CRCPoly = obj.CRC24C;
            case '16'
                CRCPoly = obj.CRC16;
            case '11'
                CRCPoly = obj.CRC11;
            case '8'
               CRCPoly = obj.CRC8;
        end
        LengthInputBit = length(InputBits);
        LengthCRCPoly = length(CRCPoly) - 1;  
        
        if InputBits(1) == 1
            XORED = [InputBits zeros(1,LengthCRCPoly)] + [CRCPoly zeros(1,LengthInputBit - 1)];
        else
            XORED = [InputBits zeros(1,LengthCRCPoly)];
        end
        for ii = 2:LengthInputBit
            if XORED(ii)*0.5 ~= floor(XORED(ii)*0.5)
                XORED(ii:ii+LengthCRCPoly) = XORED(ii:ii+LengthCRCPoly) + CRCPoly;
            end
        end
        CRCAttached = [InputBits mod(XORED(end-LengthCRCPoly+1:end),2)];        
    end
%%
    function CRCDetected = Detect(varargin)
        % Performs CRC error detection calculation and returns the result
        obj = Coding.ChannelCodingComponents.CRC;
        InputBits = varargin{1};
        switch varargin{2}
            case '24A'
                CRCPoly = obj.CRC24A;
            case '24B'
                CRCPoly = obj.CRC24B;
            case '24C'
                CRCPoly = obj.CRC24C;
            case '16'
                CRCPoly = obj.CRC16;
            case '11'
                CRCPoly = obj.CRC11;
            case '8'
               CRCPoly = obj.CRC8;
        end
        LengthCRCPoly = length(CRCPoly) - 1;  
        LengthInputBit = length(varargin{1})-LengthCRCPoly;

        if InputBits(1) == 1
            XORED = InputBits + [CRCPoly zeros(1,LengthInputBit - 1)];
        else
            XORED = InputBits;
        end
        for ii = 2:LengthInputBit
            if XORED(ii)*0.5 ~= floor(XORED(ii)*0.5)
                XORED(ii:ii+LengthCRCPoly) = XORED(ii:ii+LengthCRCPoly) + CRCPoly;
            end
        end
        CRCDetected = sum(mod(XORED(end-LengthCRCPoly+1:end),2)) > 0;        
    end
end
end