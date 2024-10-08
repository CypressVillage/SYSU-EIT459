classdef CircularBuffer < handle
% CircularBuffer Class (CircularBuffer)
%  
%   Author
%       - Bashar Tahir, btahir@nt.tuwien.ac.at
%         (c) 2018 Institute of Telecommunications, TU Wien.
%         www.nt.tuwien.ac.at 
% 
%   This class performs Circular Buffer Rate Matching (CBRM). Furthermore,
%   it supports Limited Buffer Rate Matching (LBRM).
%   Chase Combining (CC) and Incremental Redundancy (IR) are supported.
%   
%   HARQ can be disabled by settting ChannelCoding.HARQType to 'CC' and
%   making sure the ChannelCoding.RV = 0

%% LDPC related properties
properties (Constant)
    BGk = [0, 17/66, 33/66, 56/66
           0, 13/50, 25/50, 43/50];
end

%%
methods (Static)
    
    function outSeq = generateSequence(inputSeq, Block, CodingScheme, LDPCBaseGraph, LDPCZc, outputLength, bufferLength, SoftBufferRatio, subBlockInterleaverRowLength, HARQType, RV)
        outSeq = zeros(1, outputLength); 
        NULL = 98811220033.445911;
        
        Block.NullsMap = find(inputSeq==NULL);
        softBufferLength = ceil(SoftBufferRatio*bufferLength);
        R = subBlockInterleaverRowLength;
        switch HARQType
            case 'CC'
                switch CodingScheme
                    case 'TB-Convolutional'
                        k0 = 0;
                    case 'Turbo'
                        k0 = 2*R;
                    case 'LDPC'
                        k0 = 0;
                    case 'Polar'
                        k0 = 0;
                end      
            case 'IR'
                % Only Turbo and LDPC support IR
                switch CodingScheme
                    case 'Turbo'
                        k0 = R*(2*ceil(softBufferLength/(8*R))*RV + 2);
                    case 'LDPC'
                        k0 = floor(obj.BGk(LDPCBaseGraph, RV+1)*softBufferLength/LDPCZc)*LDPCZc;
                    otherwise
                        error('Only Turbo and LDPC coding support Incremental Redundancy (IR) HARQ.');
                end     
            otherwise
                error('The HARQ type is not valid.');
        end 
        Block.k0 = k0;
        Block.SoftBufferLength = softBufferLength;
        
        % Bits selection
        j = 0;
        k = 1;
        while k <= outputLength
            ind = mod(j+k0,softBufferLength) + 1;
            if inputSeq(ind) ~= NULL
                outSeq(k) = inputSeq(ind);
                k = k + 1;
            end
            j = j + 1 ;
        end
    end
    
    
    function CombinedLLRs = CombineLLRs(d, Block)
        % Performs LLRs combining and H-ARQ operation (TO DO)
        
        NULL = 98811220033.445911;
        
%         if isempty(Block.HARQBuffer) 
            % Nulls reinsertion 
            CombinedLLRs = zeros(1, Block.BaseBufferLength);
            CombinedLLRs(Block.NullsMap) = NULL;
%         else
%             CombinedLLRs = Block.HARQBuffer;
%         end
        
        j = 1;
        k = Block.k0;
        while j <= Block.MatchedLength 
            ind = mod(k,Block.SoftBufferLength)+1;
            if CombinedLLRs(ind) ~= NULL
                CombinedLLRs(ind) = CombinedLLRs(ind) + d(j);
                j = j + 1;
            end
            k = k + 1;
        end
%         Block.HARQBuffer = CombinedLLRs;
    end
end    
end

