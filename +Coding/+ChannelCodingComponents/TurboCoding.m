classdef TurboCoding < handle
% Turbo Coding Object (TurboCoding)
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
%       - Maximum A Posteriori (MAP) iterative decoder
%          -> Log-MAP
%          -> MAX-Log-MAP
%          -> Linear-Log-MAP

properties
    InputBlockSizes      % Supported code block sizes
    NrCodeBlocks         % Number of code blocks
    codeRate             % Code rate
    CodedBitsLength      % Total length of coded bits at the encoder output
    RateMatching         % The Rate Matching object
    BaseCodeRate = 1/3;  % Base code rate of the this coding scheme
    DecodingAlgorithm    % The Decoding Algorithm
    Iterations           % Number of decoding iterations
end

% Class specific properties
properties (SetAccess = private, Hidden = false)
    
    % Tables for Turbo Interleaving/Deinterleaving
    TurboInterF1 = [3,7,19,7,7,11,5,11,7,41,103,15,9,17,9,21,101,21,57,23,13,27,11,27,85,29,33,15,17,33,103,19,19,37,19,21,21,115,193,21,133,81,45,23,243,151,155,25,51,47,91,29,29,247,29,89,91,157,55,31,17,35,227,65,19,37,41,39,185,43,21,155,79,139,23,217,25,17,127,25,239,17,137,215,29,15,147,29,59,65,55,31,17,171,67,35,19,39,19,199,21,211,21,43,149,45,49,71,13,17,25,183,55,127,27,29,29,57,45,31,59,185,113,31,17,171,209,253,367,265,181,39,27,127,143,43,29,45,157,47,13,111,443,51,51,451,257,57,313,271,179,331,363,375,127,31,33,43,33,477,35,233,357,337,37,71,71,37,39,127,39,39,31,113,41,251,43,21,43,45,45,161,89,323,47,23,47,263];    % Turbo QPP (De)Interleaver f1 values
    TurboInterF2 = [10,12,42,16,18,20,22,24,26,84,90,32,34,108,38,120,84,44,46,48,50,52,36,56,58,60,62,32,198,68,210,36,74,76,78,120,82,84,86,44,90,46,94,48,98,40,102,52,106,72,110,168,114,58,118,180,122,62,84,64,66,68,420,96,74,76,234,80,82,252,86,44,120,92,94,48,98,80,102,52,106,48,110,112,114,58,118,60,122,124,84,64,66,204,140,72,74,76,78,240,82,252,86,88,60,92,846,48,28,80,102,104,954,96,110,112,114,116,354,120,610,124,420,64,66,136,420,216,444,456,468,80,164,504,172,88,300,92,188,96,28,240,204,104,212,192,220,336,228,232,236,120,244,248,168,64,130,264,134,408,138,280,142,480,146,444,120,152,462,234,158,80,96,902,166,336,170,86,174,176,178,120,182,184,186,94,190,480];    % Turbo QPP (De)Interleaver f2 values
   
    % States transitions and corresponding output (used by the decoder)
    
    TurboCodingCurrStates = [0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7];      % Current states of the State Diagram.
    TurboCodingNextStates = [0,4,4,0,5,1,1,5,2,6,6,2,7,3,3,7];      % Next states of the State Diagram.
    TurboCodingTransitionOut = [1,1,-1,-1,1,1,-1,-1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1];      % Ouput of the corresponding state transition
end
 
methods
%%
    function obj = TurboCoding(varargin)
    % Construct a new TurboCoding object
         if numel(varargin) == 1
              obj.NrCodeBlocks = varargin{1}.NrCodeBlocks;
              obj.DecodingAlgorithm = varargin{1}.DecodingAlgorithm;
              obj.codeRate = varargin{1}.codeRate;
              obj.CodedBitsLength = varargin{1}.CodedBitsLength;
              obj.RateMatching = varargin{1}.RateMatching;
              obj.InputBlockSizes = varargin{1}.InputBlockSizes;
              obj.Iterations = varargin{1}.Iterations;
          else
              error('Number of input variables must be 1');
          end
    end
%%
    function CodedBitsBlock = Encode(varargin)
        % Turbo encode the input bits block
        %
        % CodedBits = Turbo.Encode(Input bits)
        obj = varargin{1};
        BlockBits = varargin{2};
        Block = varargin{3};
        
        % Parameters
        K = Block.InputBlockSize;
        D = K + 4;  % Four tailing bits
        % Initilize first encoder
        e1_s0 = 0;
        e1_s1 = 0;
        e1_s2 = 0;
        % Initilize second encoder
        e2_s0 = 0;
        e2_s1 = 0;
        e2_s2 = 0;
        % Streams
        d0 = zeros(1,D);
        d1 = zeros(1,D);
        d2 = zeros(1,D);
        xk = zeros(1,3);
        zk = zeros(1,3);
        zdk = zeros(1,3);
        xdk = zeros(1,3);
        % Turbo Interleaving
        BlockBitsInter = obj.TurboInterleaver(BlockBits, find(obj.InputBlockSizes == K));
        % Coded streams generation
        inputS = 0;
        for k = 1:K
            e1_ci = BlockBits(k);
            e2_ci = BlockBitsInter(k);
            d0(k) = e1_ci;

            inputS = xor(e1_ci,xor(e1_s1,e1_s2));
            d1(k) = xor(e1_s2,xor(e1_s0,inputS));   
            e1_s2 = e1_s1;
            e1_s1 = e1_s0;
            e1_s0 = inputS;

            inputS = xor(e2_ci,xor(e2_s1,e2_s2));
            d2(k) = xor(e2_s2,xor(e2_s0,inputS)); 
            e2_s2 = e2_s1;  
            e2_s1 = e2_s0;
            e2_s0 = inputS;
        end
        % Trellis Termination
        for k = 1:3
            e1_ci = xor(e1_s1,e1_s2);
            e2_ci = xor(e2_s1,e2_s2);      
            xk(k) = e1_ci;
            zk(k) = xor(e1_s2,xor(e1_s0,xor(e1_ci,xor(e1_s1,e1_s2))));
            zdk(k) = xor(e2_s2,xor(e2_s0,xor(e2_ci,xor(e2_s1,e2_s2))));      
            xdk(k) = e2_ci;      
            e1_s2 = e1_s1;
            e1_s1 = e1_s0;
            e1_s0 = 0;       
            e2_s2 = e2_s1;  
            e2_s1 = e2_s0;
            e2_s0 = 0;
        end
        d0(K+1) = xk(1);
        d0(K+2) = zk(2);
        d0(K+3) = xdk(1);
        d0(K+4) = zdk(2);
        d1(K+1) = zk(1);
        d1(K+2) = xk(3);
        d1(K+3) = zdk(1);
        d1(K+4) = xdk(3);  
        d2(K+1) = xk(2);
        d2(K+2) = zk(3);
        d2(K+3) = xdk(2);
        d2(K+4) = zdk(3);
        
        % Take care of filler bits
        if Block.NrFillerBits > 0
            NULL = 98811220033.445911;
            d0(1:Block.NrFillerBits) = NULL;
            d1(1:Block.NrFillerBits) = NULL;
        end
        
        % Apply Rate matching
        CodedBitsBlock = obj.RateMatching.RateMatch('Turbo', d0, d1, d2, Block);
    end
 %%   
    function DecodedBlockBits = Decode(varargin)
    % Turbo decode the input Log-Likelihood Ratios (LLRs) block using the specified Decoding Algorithm
    % Apply the specified Decoding Algorithm and iteratively decode each
    % block. After each block is decoded, CRC is performed and the error
    % detection result is stored in the corresponding Block class.
    %
    % DecodedBits = Turbo.Decode(Channel LLRs)
        obj = varargin{1};  
        Block = varargin{3}; 
        
        % Apply Rate Dematching
        RateDematchedLLRs = obj.RateMatching.RateDematch(varargin{2}, Block)./2;

        d0 = RateDematchedLLRs(1:end/3);
        d1 = RateDematchedLLRs((end/3 + 1):2*end/3);
        d2 = RateDematchedLLRs((2*end/3 + 1):end);

        % Iterative decoding
        TurboInterleaverIndex = find(obj.InputBlockSizes == Block.InputBlockSize);
        d0tInter = obj.TurboInterleaver(d0(1:end-4), TurboInterleaverIndex);
        L21_inter = zeros(1,length(d0)-4);

        % Trellis Termination

        xk(1) = d0(end-3);
        xk(2) = d2(end-3);
        xk(3) = d1(end-2);

        zk(1) = d1(end-3);
        zk(2) = d0(end-2);
        zk(3) = d2(end-2);

        zkd(1) = d1(end-1);
        zkd(2) = d0(end);
        zkd(3) = d2(end);

        xkd(1) = d0(end-1);
        xkd(2) = d2(end-1);
        xkd(3) = d1(end);

        d0 = [d0(1:end-4) xk];
        d1 = [d1(1:end-4) zk];
        d2 = [d2(1:end-4) zkd];
        d0tInter = [d0tInter xkd];
        
        for i = 1:obj.Iterations
            L12 = obj.SISODecode(d0, d1, [L21_inter 0 0 0], 1);
            L12_inter = obj.TurboInterleaver(L12(1:end-3), TurboInterleaverIndex);
            L21 = obj.SISODecode(d0tInter, d2, [L12_inter 0 0 0], 2);
            L21_inter = obj.TurboDeinterleaver(L21(1:end-3), TurboInterleaverIndex);
        end
        TotalLLRs = 2.*d0 + [L21_inter 0 0 0] + L12;

        % Bits decision
        DecodedBlockBits = (1-sign(TotalLLRs(1:end-3)))./2;

        % Subblock CRC calculation
        if obj.NrCodeBlocks  ~= 1
            Block.CRCDetectionResult = Coding.ChannelCodingComponents.CRC.Detect(DecodedBlockBits, '24B');
            DecodedBlockBits = DecodedBlockBits(1:end-24);
        end 
    end
 %%
    function outputLLRs = SISODecode(varargin)
    % Performs the Soft-Input Soft-Ouput (BCJR) Log-MAP decoding 
    % 
    % outputLLRs = TurboCoding.SISODecode(1st stream, 2nd stream, extrinsic Stream, Decoder ID) 
        obj = varargin{1};
  
    % BCJR SISO decoding
        switch obj.DecodingAlgorithm
            case 'Log-MAP'
                outputLLRs = Coding.ChannelCodingComponents.turboDecodeMEX(varargin{2}, varargin{3}, varargin{4}, varargin{5}, obj.TurboCodingCurrStates, obj.TurboCodingNextStates, obj.TurboCodingTransitionOut, 0);
            case 'MAX-Log-MAP'
                outputLLRs = Coding.ChannelCodingComponents.turboDecodeMEX(varargin{2}, varargin{3}, varargin{4}, varargin{5}, obj.TurboCodingCurrStates, obj.TurboCodingNextStates, obj.TurboCodingTransitionOut, 1);
            case 'Linear-Log-MAP'
                outputLLRs = Coding.ChannelCodingComponents.turboDecodeMEX(varargin{2}, varargin{3}, varargin{4}, varargin{5}, obj.TurboCodingCurrStates, obj.TurboCodingNextStates, obj.TurboCodingTransitionOut, 2);
            otherwise
                error('The decoding algorithm is not supported.');
        end  
    end   
end
   

methods
%%  
    function InterleavedBits = TurboInterleaver(varargin)
        % QPP Turbo Interleaver
        obj = varargin{1};
        Index = varargin{3};
        Input = varargin{2};
        K = obj.InputBlockSizes(Index);
        InterleavedBits = zeros(1, length(Input));
        f1 = obj.TurboInterF1(Index);
        f2 = obj.TurboInterF2(Index);
        
        for k = 0:K-1
            InterleavedBits(k+1) = Input(mod(f1*k+f2*k^2,K)+1);
        end
    end
%%
    function InterleavedBits = TurboDeinterleaver(varargin)
        % QPP Turbo Deinterleaver
        obj = varargin{1};
        Index = varargin{3};
        Input = varargin{2};
        K = obj.InputBlockSizes(Index);
        InterleavedBits = zeros(1, length(Input));
        f1 = obj.TurboInterF1(Index);
        f2 = obj.TurboInterF2(Index);
        for k = 0:K-1
            InterleavedBits(mod(f1*k+f2*k^2,K)+1) = Input(k+1);
        end
    end   
end
end

