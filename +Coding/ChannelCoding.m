classdef ChannelCoding < handle
% Channel Coding Object (ChannelCoding)
% 
%   Author
%       - Bashar Tahir, btahir@nt.tuwien.ac.at
%         (c) 2018 Institute of Telecommunications, TU Wien.
%         www.nt.tuwien.ac.at 
%  
%   The object supports the following coding schemes
%       - Tail-Baiting Convolutional Code, LTE-A encoding.
%       - Turbo Code, LTE-A encoding.
%       - Low Density Parity Check (LDPC) Code, 5G-NR encoding.
%       - Polar Code, custom construction.
%     
%   Decoding Algorithms
%       - TB Convolutional and Turbo codes: Maximum A Posteriori (Log-MAP, MAX-Log-MAP, Linear-Log-MAP).
%       - LDPC code: Layered Belief Propagation (Sum-Product, Min-Sum, PWL-Min-Sum).
%       - Polar code: Successive Cancellation (SC, List-SC, CRC-List-SC).

%% Public properties
properties (SetAccess = public)
    CodingScheme         % The Coding Scheme: 'TB-Convolutional', 'Turbo', 'LDPC', 'Polar'
    DecodingAlgorithm    % The Decoding Algorithm: 'Log-MAP', 'MAX-Log-MAP', 'Linear-Log-MAP', 'Sum-Product', 'Min-Sum', 'PWL-Min-Sum', 'SC', 'List-SC', 'CRC-List-SC'
    codeRate             % Code rate
    CodedBitsLength      % Total length of coded bits at the encoder output
    CRCPoly              % Outter CRC polynomial
    CRCLength            % Outter CRC length
    CRCDetectionResult   % Outter CRC detection result
    HARQType = 'CC'      % HAQR Type: Chase Combining (CC) or Incremental Redundancy (IR)
    RV = 0;              % Current Redundancy Version
    SoftBufferRatio      % Control the length of rate matching soft buffer.
    ModulationOrder      % The used modulation order.
end

%% Read-only public properties
properties (SetAccess = private, Hidden = false)
    NrCodeBlocks         % Number of code blocks
    Blocks = [];         % Blocks Container based on the Block class
    BaseCodeRate = 1     % Base code rate
    Iterations           % Number of decoding iterations   
end

%% CodingComponents objects
properties (SetAccess = private, Hidden = false)
    TBConvCoding         % The TBConvCoding object
    TurboCoding          % The TurboCoding object
    LDPCCoding           % The LDPCCoding object
    PolarCoding          % The PolarCoding object  
    RateMatching         % The RateMatching object 
end

%% Turbo and Convolutional Code properties
properties
    InputBlockSizes = [40,48,56,64,72,80,88,96,104,112,120,128,136,144,152,160,168,176,184,192,200,208,216,224,232,240,248,256,264,272,280,288,296,304,312,320,328,336,344,352,360,368,376,384,392,400,408,416,424,432,440,448,456,464,472,480,488,496,504,512,528,544,560,576,592,608,624,640,656,672,688,704,720,736,752,768,784,800,816,832,848,864,880,896,912,928,944,960,976,992,1008,1024,1056,1088,1120,1152,1184,1216,1248,1280,1312,1344,1376,1408,1440,1472,1504,1536,1568,1600,1632,1664,1696,1728,1760,1792,1824,1856,1888,1920,1952,1984,2016,2048,2112,2176,2240,2304,2368,2432,2496,2560,2624,2688,2752,2816,2880,2944,3008,3072,3136,3200,3264,3328,3392,3456,3520,3584,3648,3712,3776,3840,3904,3968,4032,4096,4160,4224,4288,4352,4416,4480,4544,4608,4672,4736,4800,4864,4928,4992,5056,5120,5184,5248,5312,5376,5440,5504,5568,5632,5696,5760,5824,5888,5952,6016,6080,6144];   % Supported input codeblock sizes (Based on LTE-A)   
end

%% LDPC Code properties
properties
    LDPCBaseGraph        % Indicates which LDPC basegraph is used
    LDPCZc               % Indicates which extension factor is used.
end

%% Polar Code properties
properties
    F1024                 % Kronecker Matrix for Polar code of length n = 1024
    F2048                 % Kronecker Matrix for Polar code of length n = 2048
    F4096                 % Kronecker Matrix for Polar code of length n = 4096
    F8192                 % Kronecker Matrix for Polar code of length n = 8192
    F16384                % Kronecker Matrix for Polar code of length n = 16384
    PolarDefaultBaseCodeRate = 1/3;     % Polar base code rate
    ListSize              % The list size for the polar list decoder
end

%% Public Methods
methods (Access = public, Hidden = false) 
%%
    function obj = ChannelCoding(varargin)
    % Construct a new ChannelCoding object
    % ChannelCoding1 = ChannelCoding(CodingScheme, DecodingAlgorithm,
    % CodeRate, Iterations)
            obj.CodingScheme = varargin{1};
            obj.DecodingAlgorithm = varargin{2};
            obj.codeRate = varargin{3}; 
            switch obj.CodingScheme
                case 'TB-Convolutional'
                    obj.BaseCodeRate = 1/3;
                case 'Turbo'
                    obj.BaseCodeRate = 1/3;
                    obj.Iterations = varargin{4};
                case 'LDPC'
%                     obj.BaseCodeRate = obj.LDPCDefaultBaseCodeRate;
                    obj.Iterations = varargin{4};
                case 'Polar'
                    obj.BaseCodeRate = obj.codeRate;
                    load('Coding/ChannelCodingComponents/PolarPrecalculated.mat');
                    obj.F1024 = F1024;
                    obj.F2048 = F2048;
                    obj.F4096 = F4096;
                    obj.F8192 = F8192;
                    obj.F16384 = F16384;
                    clearvars F1024 F2048 F4096 F8192 F16384
                    obj.InputBlockSizes =  obj.InputBlockSizes(1:177); 
                    if strcmp(obj.DecodingAlgorithm, 'SC') == 0
                        obj.ListSize = varargin{4};
                    end
                otherwise
                    error('The coding scheme is not valid.');
            end 
            
    end
%%
    function CodedBits = encode(varargin)
    % encode Encode the input bits
    % Pass the input bits to the specific encoder based on the selected
    % CodingScheme, and updates the HARQ parameters.
    %
    % CodedBits = ChannelObject1.Decode(Input bits)
        obj = varargin{1};
        obj.RateMatching.HARQType = obj.HARQType;
        
        InputBits = varargin{2}';
        InputBits = Coding.ChannelCodingComponents.CRC.Attach(InputBits, obj.CRCPoly);

        CodedBits = [];
        BlockStartIndex = 1;
        
        
        % Add filler bits if needed for TB-Convolutional and Turbo.
        switch obj.CodingScheme
            case {'TB-Convolutional', 'Turbo', 'Polar'}
                InputBits = [zeros(1, obj.Blocks(1).NrFillerBits) InputBits];
        end
 
        for CBi = 1:obj.NrCodeBlocks
            switch obj.CodingScheme
                case {'TB-Convolutional', 'Turbo', 'Polar'}
                    if obj.NrCodeBlocks == 1
                        BlockBits = InputBits(BlockStartIndex:BlockStartIndex + obj.Blocks(CBi).InputBlockSize - 1);
                        BlockStartIndex = BlockStartIndex + obj.Blocks(CBi).InputBlockSize;
                    else
                        BlockBits = Coding.ChannelCodingComponents.CRC.Attach(InputBits(BlockStartIndex:BlockStartIndex + obj.Blocks(CBi).InputBlockSize - 25), '24B');
                        BlockStartIndex = BlockStartIndex + obj.Blocks(CBi).InputBlockSize - 24;
                    end
                case 'LDPC'
                    if obj.NrCodeBlocks == 1
                        BlockBits = [InputBits(BlockStartIndex:BlockStartIndex + obj.Blocks(CBi).InputBlockSize - obj.Blocks(CBi).NrFillerBits - 1) zeros(1,obj.Blocks(CBi).NrFillerBits)];
                        BlockStartIndex = BlockStartIndex + obj.Blocks(CBi).InputBlockSize;
                    else
                        BlockBits = [Coding.ChannelCodingComponents.CRC.Attach(InputBits(BlockStartIndex:BlockStartIndex + obj.Blocks(CBi).InputBlockSize - obj.Blocks(CBi).NrFillerBits - 25), '24B') zeros(1,obj.Blocks(CBi).NrFillerBits)];
                        BlockStartIndex = BlockStartIndex + obj.Blocks(CBi).InputBlockSize - obj.Blocks(CBi).NrFillerBits - 24;
                    end
            end
            
            switch obj.CodingScheme
                case 'TB-Convolutional'
                    CodedBitsBlock = obj.TBConvCoding.Encode(BlockBits, obj.Blocks(CBi));            
                case 'Turbo'
                    CodedBitsBlock = obj.TurboCoding.Encode(BlockBits, obj.Blocks(CBi)); 
                case 'LDPC'
                    CodedBitsBlock = obj.LDPCCoding.Encode(BlockBits, obj.Blocks(CBi)); 
                case 'Polar'
                    CodedBitsBlock = obj.PolarCoding.Encode(BlockBits, obj.Blocks(CBi)); 
                otherwise
                    error('The coding scheme is not valid.');
            end
            
            % Concatenate the coded blocks
            CodedBits = [CodedBits CodedBitsBlock];
        end
        CodedBits = transpose(CodedBits);
    end
 %%   
    function DecodedBits = decode(varargin)
    % decode Decode the input Log-Likelihood Ratios (LLRs)
    % Passes the input LLRs to the specific decoder based on the selected
    % CodingScheme, and performs outter CRC error detection.
    %
    % DecodedBits = ChannelObject1.Decode(Channel LLRs)
        obj = varargin{1};
        UnpuncturedLength = 0;
        for CBi = 1:obj.NrCodeBlocks
            UnpuncturedLength = UnpuncturedLength + obj.Blocks(CBi).MatchedLength;
        end
        ChannelLLRs = zeros(1, UnpuncturedLength);
        ChannelLLRs(1:length(varargin{2})) = varargin{2}';
        
        DecodedBits = [];
        BlockStartIndex = 1;
    
        for CBi = 1:obj.NrCodeBlocks
            ChannelLLRsBlock = ChannelLLRs(BlockStartIndex:BlockStartIndex + obj.Blocks(CBi).MatchedLength - 1);
            BlockStartIndex = BlockStartIndex + obj.Blocks(CBi).MatchedLength;
            
            switch obj.CodingScheme
                case 'TB-Convolutional'
                    DecodedBlockBits = obj.TBConvCoding.Decode(ChannelLLRsBlock, obj.Blocks(CBi));    
                case 'Turbo'
                    DecodedBlockBits = obj.TurboCoding.Decode(ChannelLLRsBlock, obj.Blocks(CBi));    
                case 'LDPC'
                    DecodedBlockBits = obj.LDPCCoding.Decode(ChannelLLRsBlock, obj.Blocks(CBi));    
                case 'Polar'
                    DecodedBlockBits = obj.PolarCoding.Decode(ChannelLLRsBlock, obj.Blocks(CBi)); 
                otherwise
                    error('The coding scheme is not valid.');
            end 
            
            DecodedBits = [DecodedBits DecodedBlockBits];
        end

        obj.CRCDetectionResult = Coding.ChannelCodingComponents.CRC.Detect(DecodedBits, obj.CRCPoly);  
        if obj.NrCodeBlocks == 1
            obj.Blocks(1).CRCDetectionResult = obj.CRCDetectionResult;
        end
        switch obj.CodingScheme
            case {'TB-Convolutional', 'Turbo', 'Polar'}
                % Remove Filler bits
                DecodedBits = transpose(DecodedBits(obj.Blocks(1).NrFillerBits+1:end-obj.CRCLength));
            case 'LDPC'
                % Filler bits are already removed at prior stage
                DecodedBits = transpose(DecodedBits(1:end-obj.CRCLength));
        end
        
    end
 %%   
 function requiredInputLength = update(varargin)
        % Updates the object, performs code block segmentation and returns the required input length
        
        obj = varargin{1};  
        obj.codeRate = varargin{4};
        obj.ModulationOrder = varargin{5};
        obj.SoftBufferRatio = varargin{6};
        
        % Find the required input length if needed, depending on whether
        % the input parameter is the information length or code length.
        switch varargin{2}
            case 'Input'
                inputLength = varargin{3};
                requiredInputLength = []; 
            case 'Output'
                obj.CodedBitsLength = varargin{3};
                inputLength = floor(varargin{3}*obj.codeRate);
                if strcmp(obj.CodingScheme, 'Turbo')
                    inputLength = inputLength - 12;
                end
                requiredInputLength = inputLength;
            otherwise
                error(['The length indicator is either ' char(39) 'Input' char(39) ' or ' char(39) 'Output' char(39) '.']);
        end
        
        % Code block segmenation
        switch obj.CodingScheme
            case {'Turbo', 'TB-Convolutional'}
                % LTE
                if isempty(obj.CRCPoly)
                    obj.CRCPoly = '24A';
                    obj.CRCLength = 24;
                end
                B = inputLength + obj.CRCLength;
                
                Z = 6144;
                if B <= Z
                    L = 0;
                    C = 1;
                    Bp = B;
                else
                    L = 24;
                    C = ceil(B/(Z-L));
                    Bp = B + C*L;
                end
                obj.NrCodeBlocks = C;
                
                for iKplus = 1:length(obj.InputBlockSizes)
                    if C*obj.InputBlockSizes(iKplus) >= Bp
                        Kplus = obj.InputBlockSizes(iKplus);
                        break;
                    end
                end
                if C == 1
                    Cplus = 1;
                    Kminus = 0;
                    Cminus = 0;
                elseif C > 1
                    Kminus = obj.InputBlockSizes(iKplus - 1);
                    dK = Kplus - Kminus;
                    Cminus = floor((C*Kplus-Bp)/dK);
                    Cplus = C - Cminus;
                end
                F = Cplus*Kplus + Cminus*Kminus - Bp;
                
                obj.Blocks = [];
                BlkIndex = 1;
                obj.Blocks = [obj.Blocks Coding.ChannelCodingComponents.Block(obj, Kplus, obj.codeRate, obj.BaseCodeRate, BlkIndex, F)];
                BlkIndex = BlkIndex + 1;
                for i = 2:Cplus
                    obj.Blocks = [obj.Blocks Coding.ChannelCodingComponents.Block(obj, Kplus, obj.codeRate, obj.BaseCodeRate, BlkIndex, 0)];
                    BlkIndex = BlkIndex + 1;
                end
                for i = 1:Cminus
                    obj.Blocks = [obj.Blocks Coding.ChannelCodingComponents.Block(obj, Kminus, obj.codeRate, obj.BaseCodeRate, BlkIndex, 0)];
                    BlkIndex = BlkIndex + 1;
                end 
                
            case {'LDPC'} 
                % 5G-NR
                if isempty(obj.CRCPoly)
                    if inputLength > 3824
                        obj.CRCPoly = '24A';
                        obj.CRCLength = 24;
                    else
                        obj.CRCPoly = '16';
                        obj.CRCLength = 16;
                    end
                end
                B = inputLength + obj.CRCLength;
                
                if (inputLength <= 292) || (inputLength <=3824 && obj.codeRate <= 2/3) || (obj.codeRate <=1/4)
                    obj.LDPCBaseGraph = 2;
                else
                    obj.LDPCBaseGraph = 1;
                end
                
                if obj.LDPCBaseGraph == 1
                   obj.BaseCodeRate = 1/3;
                   Kcb = 8448;
                else
                   obj.BaseCodeRate = 1/5;
                   Kcb = 3840;
                end
                
                if B <= Kcb
                    L = 0;
                    C = 1;
                    Bp = B;
                else
                    L = 24;
                    C = ceil(B/(Kcb-L));
                    Bp = B + C*L;
                end
                obj.NrCodeBlocks = C;
                
                Kp = Bp/C;
                if obj.LDPCBaseGraph == 1
                    Kb = 22;
                else
                    if B > 640
                        Kb = 10;
                    elseif B > 560
                        Kb = 9;
                    elseif B > 192
                        Kb = 8;
                    else
                        Kb = 6;
                    end
                end
                Zlist = Coding.ChannelCodingComponents.LDPCCoding.LiftingSizeSets;
                Zc = min(Zlist(Kb.*Zlist>=Kp)); 
                obj.LDPCZc = Zc;
                [iSet,~] = find(Zlist==Zc);
                if obj.LDPCBaseGraph == 1
                    K = 22*Zc;
                else
                    K = 10*Zc;
                end
                obj.Blocks = [];
                for r = 0:C-1
                    NrFillerBits = K-Kp;
                    obj.Blocks = [obj.Blocks Coding.ChannelCodingComponents.Block(obj, Kp + NrFillerBits, obj.codeRate, obj.BaseCodeRate, r + 1, NrFillerBits)];
                end
            case {'Polar'}
                if isempty(obj.CRCPoly)
                    obj.CRCPoly = '24A';
                    obj.CRCLength = 24;
                end
                B = inputLength + obj.CRCLength;
                
                Z = 6144;
                if B <= Z
                    L = 0;
                    C = 1;
                    Bp = B;
                else
                    L = 24;
                    C = ceil(B/(Z-L));
                    Bp = B + C*L;
                end
                obj.NrCodeBlocks = C;
                
                Kplus = ceil(Bp/C);
                Cplus = C;
                
                obj.Blocks = [];
                BlkIndex = 1;
                obj.Blocks = [obj.Blocks Coding.ChannelCodingComponents.Block(obj, Kplus, obj.codeRate, obj.BaseCodeRate, BlkIndex, 0)];
                BlkIndex = BlkIndex + 1;
                for i = 2:Cplus
                    obj.Blocks = [obj.Blocks Coding.ChannelCodingComponents.Block(obj, Kplus, obj.codeRate, obj.BaseCodeRate, BlkIndex, 0)];
                    BlkIndex = BlkIndex + 1;
                end
        end
        switch varargin{2}
            case 'Input'
                obj.CodedBitsLength = 0;
                for iCB = 1:obj.NrCodeBlocks
                    obj.CodedBitsLength = obj.CodedBitsLength + obj.Blocks(iCB).CodedBlockSize;
                end
            case 'Output'
            otherwise
                error(['The length indicator is either ' char(39) 'Input' char(39) ' or ' char(39) 'Output' char(39) '.']);
        end
        obj.RateMatching = Coding.ChannelCodingComponents.RateMatching(obj);
        switch obj.CodingScheme
            case 'TB-Convolutional'       
                obj.TBConvCoding = Coding.ChannelCodingComponents.TBConvCoding(obj);
            case 'Turbo'
                obj.TurboCoding = Coding.ChannelCodingComponents.TurboCoding(obj);
            case 'LDPC'
                obj.LDPCCoding = Coding.ChannelCodingComponents.LDPCCoding(obj, iSet, obj.LDPCBaseGraph, Zc);
            case 'Polar'
                Coding.ChannelCodingComponents.PolarCoding.ConstructPolarCode(obj);
                if obj.NrCodeBlocks > 1
                    PolarCodeCRC = '24B';
                else
                    PolarCodeCRC = obj.CRCPoly;
                end
                obj.PolarCoding = Coding.ChannelCodingComponents.PolarCoding(obj, PolarCodeCRC);
            otherwise
                error('The coding scheme is not valid.');
        end
    end  
end
end
