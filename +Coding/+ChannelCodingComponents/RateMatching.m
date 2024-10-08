classdef RateMatching < handle
%                  Rate Matching Object (RateMatching)
%   Author
%       - Bashar Tahir, btahir@nt.tuwien.ac.at
%         (c) 2018 Institute of Telecommunications, TU Wien.
%         www.nt.tuwien.ac.at 
%  
%   Encoder side
%       - Rate Matcher
% 
%   Decoder side
%       - Rate Dematcher

properties
    NrCodeBlocks         % Number of code blocks
    Blocks               % Blocks Container based on the Block class
    codeRate             % Code rate
    CodedBitsLength      % Total length of coded bits at the encoder output
    BaseCodeRate = 1/3;  % Base code rate of the this coding scheme
    CodingScheme         % The used coding scheme
    LDPCBaseGraph
    LDPCZc
    CircularBuffer = Coding.ChannelCodingComponents.CircularBuffer   % The CircularBuffer object
    HARQType             % HAQR Type: Chase Combining (CC) or Incremental Redundancy (IR)
    RV                   % Redundancy Version of the transmission
    ModulationOrder      % The modulation order used
    SoftBufferRatio      % Controls the length of the soft buffer limitation
end

properties (Constant)
    NULL = 98811220033.445911;
end

properties (Constant)
    % Interleavers and Deinterleavers tables
    
    TBConvCodingInterleaverP = [2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31];
    TBConvCodingDeinterleaverP = [17,1,25,9,21,5,29,13,19,3,27,11,23,7,31,15,18,2,26,10,22,6,30,14,20,4,28,12,24,8,32,16]; 
    TurboSBInterd0d1P = [1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32];
    TurboSBInterd2P = [0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31];
    TurboSBDeinterd0d1P = [1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32];
    TurboSBDeinterd2P = [32,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16];
end

methods
%%
    function obj = RateMatching(varargin)
    % Construct a new RateMatching object
        if numel(varargin) == 1
              obj.NrCodeBlocks = varargin{1}.NrCodeBlocks;
              obj.Blocks = varargin{1}.Blocks;
              obj.codeRate = varargin{1}.codeRate;
              obj.CodedBitsLength = varargin{1}.CodedBitsLength;
              obj.BaseCodeRate = varargin{1}.BaseCodeRate;
              obj.HARQType = varargin{1}.HARQType;
              obj.ModulationOrder = varargin{1}.ModulationOrder;
              obj.SoftBufferRatio = varargin{1}.SoftBufferRatio;
              obj.RV = varargin{1}.RV;
              obj.LDPCBaseGraph = varargin{1}.LDPCBaseGraph;
              obj.LDPCZc = varargin{1}.LDPCZc;
        else
              error('Number of input variables must be 1');
        end
    end
%%
    function MatchedBits = RateMatch(varargin)
     % Performs Rate Matching
     % Performs subblock interleaving, bits collection and HARQ based bit selection.
     %
     % MatchedBits = RateMatching.RateMatch(CodingScheme, 1st Stream, 2nd ...
     % Stream, 3rd Stream, Block size, Index of current Block)
     %
     % For Tail-Biting Convolutional and Turbo coding schemes, all the
     % three streams are used. For LDPC, only the first two, and for Polar
     % only the first one.
     
        obj = varargin{1};
        obj.CodingScheme = varargin{2};
        Block = varargin{6};
        D = Block.InputBlockSize;
        if strcmp(obj.CodingScheme, 'Turbo') == 1
            D = D + 4;
        end
        K_PI = [];
        R = [];
        
        % Determine the rate matched output length for each block
        if  Block.Index ~= obj.NrCodeBlocks
            E = Block.CodedBlockSize;            
        else
            if obj.NrCodeBlocks == 1
                 E = obj.CodedBitsLength;
            else
                 sumT = 0;
                 for CBi = 1:obj.NrCodeBlocks - 1
                    sumT = sumT + obj.Blocks(CBi).CodedBlockSize;
                 end
                 E = obj.CodedBitsLength - sumT;
            end    
        end                       
        
        switch obj.CodingScheme
            case 'TB-Convolutional'
            % Sub-Block Interleaving
                d0 = obj.TBConvCodingsubBlkInterleaver(varargin{3});
                d1 = obj.TBConvCodingsubBlkInterleaver(varargin{4});
                d2 = obj.TBConvCodingsubBlkInterleaver(varargin{5});
            % Bits Collection
                C = 32;                                     
                R = ceil(D/C);                              
                K_PI = R*C;                                 
                K_w = 3*K_PI;                               
                w = [d0 d1 d2];                             
            case 'Turbo'
            % Sub-Block Interleaving
                d0 = obj.TurbosubBlkInterleaverd0d1(varargin{3});
                d1 = obj.TurbosubBlkInterleaverd0d1(varargin{4});
                d2 = obj.TurbosubBlkInterleaverd2(varargin{5});
            % Bits Collection
                C = 32;                                     
                R = ceil(D/C);                              
                K_PI = R*C;                                 
                K_w = 3*K_PI;                               
                w = zeros(1, K_w);   
                w(1:K_PI) = d0;
                w(K_PI+1:2:end) = d1;
                w(K_PI+2:2:end) = d2;
            case 'LDPC'
            % No interleaving at this stage for 5G-NR LDPC codes
                w = varargin{3};
                K_w = length(w);   
            case 'Polar'
            % Interleaving    
            if obj.codeRate >= obj.BaseCodeRate
                w = obj.PolarInterleaver(varargin{3}(1:Block.CodedBlockSize));
            else
                w = obj.PolarInterleaver(varargin{3}(1:Block.BaseCodedBlockSize));             
            end
            K_w = length(w);
            otherwise
                error('The coding scheme is not valid.');
        end
        Block.K_PI = K_PI;
        Block.BaseBufferLength = K_w;
        Block.MatchedLength = E;
    % Bit Selection
        % Circular-buffer processing
        e = obj.CircularBuffer.generateSequence(w, Block, obj.CodingScheme, obj.LDPCBaseGraph, obj.LDPCZc, E, K_w, obj.SoftBufferRatio, R, obj.HARQType, obj.RV);
        
        % Post circular-buffer interleaving
        switch obj.CodingScheme
            case 'LDPC'
                % 5G-NR Interleaveing
                e = obj.LDPCInterleaver(e, E);
        end
        
        MatchedBits = e;
     end
%%
    function Dematched = RateDematch(varargin)
    % Performs Rate Dematching
    %
    % Dematched = RateMatching.RateDematch(Block LLRs, Index of current block)
        obj = varargin{1};
        Block = varargin{3};
        BlockSize = Block.InputBlockSize;

        switch obj.CodingScheme
            case 'TB-Convolutional'
                % Get the size of a single stream
                K_PI = Block.K_PI;

                % Nulls reinsertion, LLRs combining and H-ARQ operation
                w = obj.CircularBuffer.CombineLLRs(varargin{2}, Block);
                
                % Bits Decollection
                d0 = w(1:K_PI);
                d1 = w(K_PI+1:2*K_PI);
                d2 = w(2*K_PI+1:end);
                
                % De-Interleaving and Nulls removal
                d0 = obj.TBConvCodingsubBlkDeinterleaver(d0, BlockSize);
                d1 = obj.TBConvCodingsubBlkDeinterleaver(d1, BlockSize);
                d2 = obj.TBConvCodingsubBlkDeinterleaver(d2, BlockSize);
                
                Dematched = [d0 d1 d2]; 
                   
            case 'Turbo'
                % Get the size of a single stream
                BlockSize = BlockSize + 4;
                K_PI = Block.K_PI;
                
                % Nulls reinsertion, LLRs combining and H-ARQ operation
                w = obj.CircularBuffer.CombineLLRs(varargin{2}, Block);

                % Bits decollection
                d0 = w(1:K_PI);
                d1 = w(K_PI+1:2:end);
                d2 = w(K_PI+2:2:end);

                % De-interleaving and nulls removal
                d0 = obj.TurbosubBlkDeinterleaverd0d1(d0, BlockSize, Block.NrFillerBits);
                d1 = obj.TurbosubBlkDeinterleaverd0d1(d1, BlockSize, Block.NrFillerBits);
                d2 = obj.TurbosubBlkDeinterleaverd2(d2, BlockSize);
                
                Dematched = [d0 d1 d2];         
                
            case 'LDPC'
                % De-Interleaving
                d = obj.LDPCDeinterleaver(varargin{2}(1:Block.MatchedLength));
                d = [d zeros(1, length(varargin{2}) - Block.MatchedLength)];
                
                % Nulls reinsertion, LLRs combining and H-ARQ operation
                w = obj.CircularBuffer.CombineLLRs(d, Block);
                
                % Set Filler bits to very high positive value (+Inf)
                w(Block.NullsMap) = 10^10;
                
                Dematched = w;

            case 'Polar'
                if obj.codeRate < obj.BaseCodeRate
                    BlockSize = Block.BaseCodedBlockSize;
                else
                    BlockSize = Block.CodedBlockSize;
                end
                
                % Nulls reinsertion, LLRs combining and H-ARQ operation
                w = obj.CircularBuffer.CombineLLRs(varargin{2}, Block);

                % De-Interleaving and Nulls removal
                Dematched = obj.PolarDeinterleaver(w, BlockSize); 
            otherwise
                error('The coding scheme is not valid.');
        end
    end
    
end   
 
%% Interleavers and Deinterleavers
methods (Access = private, Hidden = true)
%%   
    function InterleavedBits = PolarInterleaver(varargin)
        obj = varargin{1};
        BlockBits = varargin{2};
        D = length(BlockBits);
        C = 32;  
        R = ceil(D/C);
        indStart = R*C-D;
        if indStart ~=0
            BlockBits = [obj.NULL*ones(1,indStart) BlockBits];
        end    
        pre = reshape(BlockBits, C, R).';
        InterleavedBits = reshape(pre(:,obj.TBConvCodingInterleaverP), 1, []); 
    end
%% 
    function InterleavedBits = TBConvCodingsubBlkInterleaver(varargin)
        obj = varargin{1};
        BlockBits = varargin{2};
        D = length(BlockBits);
        C = 32;  
        R = ceil(D/C);
        indStart = R*C-D;
        if indStart ~=0
            BlockBits = [obj.NULL*ones(1,indStart) BlockBits];
        end    
        pre = reshape(BlockBits, C, R).';
        InterleavedBits = reshape(pre(:,obj.TBConvCodingInterleaverP), 1, []); 
    end
%%
    function InterleavedBits = LDPCInterleaver(varargin)
        E = varargin{3};
        Qm = varargin{1}.ModulationOrder;
        InterleavedBits = reshape(reshape(varargin{2}, E/Qm, []).', E, []).';
    end
%%
    function InterleavedBits = TurbosubBlkInterleaverd0d1(varargin)
        obj = varargin{1};
        BlockBits = varargin{2};
        D = length(BlockBits);
        C = 32;  
        R = ceil(D/C);
        indStart = R*C-D;
        if indStart ~=0
            BlockBits = [obj.NULL*ones(1,indStart) BlockBits];
        end    
        pre = reshape(BlockBits, C, R).';
        InterleavedBits = reshape(pre(:,obj.TurboSBInterd0d1P), 1, []);     
    end
%%    
    function InterleavedBits = TurbosubBlkInterleaverd2(varargin)
        obj = varargin{1};
        BlockBits = varargin{2};
        D = length(BlockBits);
        C = 32;
        R = ceil(D/C);
        P = obj.TurboSBInterd2P;
        K_PI = R*C;
        indStart = K_PI-D;
        if indStart ~=0
            BlockBits = [obj.NULL*ones(1,indStart) BlockBits];
        end    
        k = 0:K_PI - 1;
        InterleavedBits = BlockBits(mod(P(floor(k/R)+1) + C*mod(k,R) + 1, K_PI)+1);
    end
%%
    function Deinterleaved = PolarDeinterleaver(varargin)
        obj = varargin{1};
        BlockLLRs =  varargin{2};
        D = length(BlockLLRs);
        C = 32;  
        R = ceil(D/C);
        NrInterleaveNulls = R*C-varargin{3};
        pre = reshape(BlockLLRs, R, C);
        per = reshape(pre(:,obj.TBConvCodingDeinterleaverP).', 1, []);  
        Deinterleaved = per(NrInterleaveNulls+1:end);
    end
%%
    function Deinterleaved = LDPCDeinterleaver(varargin)
        Qm = varargin{1}.ModulationOrder;
        blockLength = length(varargin{2});
        
        Deinterleaved = reshape(reshape(varargin{2}.', [], blockLength/Qm).', blockLength, []).';
    end
%%
    function Deinterleaved = TBConvCodingsubBlkDeinterleaver(varargin)
        obj = varargin{1};
        BlockLLRs =  varargin{2};
        D = length(BlockLLRs);
        C = 32;  
        R = ceil(D/C);
        NrInterleaveNulls = R*C-varargin{3};
        pre = reshape(BlockLLRs, R, C);
        per = reshape(pre(:,obj.TBConvCodingDeinterleaverP).', 1, []);  
        Deinterleaved = per(NrInterleaveNulls+1:end);
    end
%%    
    function Deinterleaved = TurbosubBlkDeinterleaverd0d1(varargin)
        obj = varargin{1};
        BlockLLRs =  varargin{2};
        D = length(BlockLLRs);
        C = 32;  
        R = ceil(D/C);
        NrInterleaveNulls = R*C-varargin{3};
        
        pre = reshape(BlockLLRs, R, C);
        per = reshape(pre(:,obj.TurboSBInterd0d1P).', 1, []);  
        
        % Set the LLRs for filler bits to high positive value (+Inf)
        per(NrInterleaveNulls+1:NrInterleaveNulls+varargin{4}) = 10^10;
        
        Deinterleaved = per(NrInterleaveNulls+1:end);
    end
       
%%       
    function Deinterleaved = TurbosubBlkDeinterleaverd2(varargin)
            obj = varargin{1};
            BlockLLRs =  varargin{2};
            D = length(BlockLLRs);
            C = 32;  
            R = ceil(D/C);
            P = obj.TurboSBInterd2P;
            K_PI = R*C;
            NrInterleaveNulls = R*C-varargin{3};
            
            k = 0:K_PI - 1;
            Deinterleaved(mod(P(floor(k/R)+1) + C*mod(k,R) + 1, K_PI)+1) = BlockLLRs;
            Deinterleaved(1:NrInterleaveNulls) = []; 
        end
    end
    
end