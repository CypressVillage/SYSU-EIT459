classdef Block < handle
% Block Object (Block)
%
%   Author
%       - Bashar Tahir, btahir@nt.tuwien.ac.at
%         (c) 2018 Institute of Telecommunications, TU Wien.
%         www.nt.tuwien.ac.at 
%  
%   Contains all the necessary information about the transmitted block.
%   Makes the decoding process much faster, since the decoder can obtain all
%   the needed information from this class, instead of having to redo the
%   calculation by itself.  

%%
properties
    InputBlockSize          % Size of the input code block including subblock CRC
    CodedBlockSize          % Size of punctured output coded block
    BaseCodedBlockSize      % Size of unpunctured output coded block
    Index                   % Block index with respect to other blocks
    CRCDetectionResult      % Result of CRC calculation at the decoder
    NrFillerBits = 0        % The number of filler bits used for this block
    NullsMap                % The location of nulls in this block
    SoftBufferLength        % The soft length of the rate matching buffer
    BaseBufferLength        % The base length of rate matching buffer
    HARQBuffer = []         % Stores the previous transmission channel LLRs
    InterleaverRowLength    % For TB-C, Turbo, Polar
    k0                      % The starting position of the current transmission
    K_PI                    % Stream size for Convolutional and Turbo
    MatchedLength           % Length of block after rate matching.
end
% TB-C properties
properties
    InitialState            % Initial state of the TB-C encoder,
end
% Turbo properties
properties 
end
% LDPC properties
properties
% Decoder Part (Tanner Graph)

    LDPCcnMat               % Variable nodes connected to Check nodes
    LDPCvnMat               % Check nodes connected to Variable Node
end
% Polar properties
properties
    Power2Size               % The butter fly size of the Polar encoder
    FrozenSet                % The construced Frozen set for the Polar code
% Encoder Part

    KroneckerMatrix          % Generator matrix for the Polar code
end
%%
methods (Access = public, Hidden = false) 
    function obj = Block(varargin)
    % Construct a new Block class
        obj.InputBlockSize = varargin{2};
        obj.Index = varargin{5};
        obj.NrFillerBits = varargin{6};
        switch varargin{1}.CodingScheme
            case 'TB-Convolutional'
                obj.CodedBlockSize = ceil(obj.InputBlockSize/varargin{3});
                obj.BaseCodedBlockSize = ceil(obj.InputBlockSize/varargin{4});
            case 'Turbo'
                obj.CodedBlockSize = ceil((obj.InputBlockSize)/varargin{3}) - 2*obj.NrFillerBits;
                obj.BaseCodedBlockSize = ceil((obj.InputBlockSize)/varargin{4}) - 2*obj.NrFillerBits;
            case 'LDPC'
                obj.CodedBlockSize = ceil(obj.InputBlockSize/varargin{3}) - obj.NrFillerBits;
                obj.BaseCodedBlockSize = ceil((obj.InputBlockSize)/varargin{4} - obj.NrFillerBits);
            case 'Polar'
                obj.CodedBlockSize = ceil(obj.InputBlockSize/varargin{3});
                obj.BaseCodedBlockSize = ceil(obj.InputBlockSize/varargin{4});
        end
                
    end
end

end

