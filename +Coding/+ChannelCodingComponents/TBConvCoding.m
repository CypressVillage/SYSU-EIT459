classdef TBConvCoding < handle
% Tail-Biting Convolutional Coding Object (TBConvCoding)
% 
%   Author
%       - Bashar Tahir, btahir@nt.tuwien.ac.at
%         (c) 2018 Institute of Telecommunications, TU Wien.
%         www.nt.tuwien.ac.at 
%  
%   Encoding
%       - LTE-A based encoder.
% 
%   Decoding
%       - Maximum A Posteriori (MAP) decoder
%          -> Log-MAP
%          -> MAX-Log-MAP

properties
    NrCodeBlocks         % Number of code blocks
    codeRate             % Code rate
    CodedBitsLength      % Total length of coded bits at the encoder output
    RateMatching         % The Rate Matching object
    BaseCodeRate = 1/3;  % Base code rate of the this coding scheme
    DecodingAlgorithm    % The Decoding Algorithm
end

% Class specific properties
properties (SetAccess = private, Hidden = false)
    % States transitions and corresponding output (used by the decoder)
    
    TBConvCodingCurrStates = [0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,41,41,42,42,43,43,44,44,45,45,46,46,47,47,48,48,49,49,50,50,51,51,52,52,53,53,54,54,55,55,56,56,57,57,58,58,59,59,60,60,61,61,62,62,63,63];      % Current states of the State Diagram.
    TBConvCodingNextStates = [0,32,0,32,1,33,1,33,2,34,2,34,3,35,3,35,4,36,4,36,5,37,5,37,6,38,6,38,7,39,7,39,8,40,8,40,9,41,9,41,10,42,10,42,11,43,11,43,12,44,12,44,13,45,13,45,14,46,14,46,15,47,15,47,16,48,16,48,17,49,17,49,18,50,18,50,19,51,19,51,20,52,20,52,21,53,21,53,22,54,22,54,23,55,23,55,24,56,24,56,25,57,25,57,26,58,26,58,27,59,27,59,28,60,28,60,29,61,29,61,30,62,30,62,31,63,31,63];      % Next states of the State Diagram.
    TBConvCodingTransitionOut = [1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,-1,-1,1,1,1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1,-1,-1,-1,1,-1,-1,1,1,1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,-1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1];      % Ouput of the corresponding state transition  
end
 
methods
%%
    function obj = TBConvCoding(varargin)
    % Construct a new TBConvCoding object
         if numel(varargin) == 1
              obj.NrCodeBlocks = varargin{1}.NrCodeBlocks;
              obj.codeRate = varargin{1}.codeRate;
              obj.CodedBitsLength = varargin{1}.CodedBitsLength;
              obj.RateMatching = varargin{1}.RateMatching;
              obj.DecodingAlgorithm = varargin{1}.DecodingAlgorithm;
          else
              error('Number of input variables must be 1');
          end
    end
%%
    function CodedBitsBlock = Encode(varargin)
        % Convolutionally encode the input bits block
        %
        % CodedBits = TBConvCoding.Encode(Input bits)

        % Store the encoder initial state in the corresponding Block class
        obj = varargin{1};
        BlockBits = varargin{2};
        Block = varargin{3};
        
        InitialState = bi2de(BlockBits(end-5:end));          
        Block.InitialState = InitialState;   
        % Parameters
        D = Block.InputBlockSize;              
        % Streams
        d0 = zeros(1,D);
        d1 = zeros(1,D);
        d2 = zeros(1,D);
        % Tail Biting
        s0 = BlockBits(end);
        s1 = BlockBits(end-1);
        s2 = BlockBits(end-2);
        s3 = BlockBits(end-3);
        s4 = BlockBits(end-4);
        s5 = BlockBits(end-5);
        ci = 0;

        % Coded streams generation   
        for k= 1:D
           ci = BlockBits(k);   
           d0(k) = mod(ci+s1+s2+s4+s5,2);
           d1(k) = mod(ci+s0+s1+s2+s5,2);
           d2(k) = mod(ci+s0+s1+s3+s5,2);
           s5 = s4;
           s4 = s3;
           s3 = s2;
           s2 = s1;
           s1 = s0;
           s0 = ci;
        end
        % Apply Rate Matching
        CodedBitsBlock = obj.RateMatching.RateMatch('TB-Convolutional', d0, d1, d2, Block); 
    end
 %%   
    function DecodedBlockBits = Decode(varargin)
    % Convolutionally decode the input Log-Likelihood Ratios (LLRs) block using the specified Decoding Algorithm
    % Apply the specified Decoding Algorithm. After each block is decoded, CRC is performed and the error
    % detection result is stored in the corresponding Block class.
    %
    % DecodedBits = TBConvCoding.Decode(Channel LLRs)
        obj = varargin{1};  
        Block = varargin{3}; 
         
        % Apply Rate Dematching
        RateDematchedLLRs = obj.RateMatching.RateDematch(varargin{2}, Block)./2;

        % SISO (BCJR) decoding
        DecodedBlockBits = obj.TBConvCodingDecode(RateDematchedLLRs, Block.InitialState); 

        % Subblock CRC calculation
        if obj.NrCodeBlocks  ~= 1
            Block.CRCDetectionResult = Coding.ChannelCodingComponents.CRC.Detect(DecodedBlockBits, '24B');
            DecodedBlockBits = DecodedBlockBits(1:end-24);
        end 
        
    end
%%
  function DecodedBits = TBConvCodingDecode(varargin)
    % Performs the BCJR Log-MAP decoding 
    %
    % decodedBits = TBConvCoding.TBConvCodingDecode(Dematched LLRs, Initial state)
        obj = varargin{1};

        switch obj.DecodingAlgorithm
          case 'Log-MAP'
              outputLLRs = Coding.ChannelCodingComponents.TBConvDecodeMEX(varargin{2}, varargin{3}, obj.TBConvCodingCurrStates, obj.TBConvCodingNextStates, obj.TBConvCodingTransitionOut, 0);
          case 'MAX-Log-MAP'
              outputLLRs = Coding.ChannelCodingComponents.TBConvDecodeMEX(varargin{2}, varargin{3}, obj.TBConvCodingCurrStates, obj.TBConvCodingNextStates, obj.TBConvCodingTransitionOut, 1);
           otherwise
                error('The decoding algorithm is not supported.');
        end
        
        DecodedBits = (1 - sign(outputLLRs))./2;
    end   
end

end

