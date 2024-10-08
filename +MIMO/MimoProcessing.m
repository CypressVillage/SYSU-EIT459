classdef MimoProcessing < handle
    % this class contains methods and properties for MIMO modulation and
    % demodulation
    
    properties
        method      % 'spatial multiplexing'
        nStreams    % number of spatial streams
        codebook    % codebook used for MIMO   
        precodingMatrix % precoding matrix defined in the scenario file
        precodingMatrixIndex % index of the precoding matrix used in a codebook based transmission 
    end
    
    methods
        function obj = MimoProcessing( mimoMethod, nSpatialStreams )
            % this is the class constructor
            obj.method      = mimoMethod;
            
            switch mimoMethod
                case 'single transmit antenna'
                    if( nSpatialStreams ~= 1)
                        error('Only a single spatial stream is supported for transmission mode ''single transmit antenna''!');
                    end
                case {'CLSM','OLSM','TxD','custom'}
                otherwise
                    error('Transmission mode unkown!');
            end
        end
        
        function precodedDataSymbols = precoding( obj, dataSymbols, nAntennas )
            % MIMO precoding of data symbols
            
            % check inputs
            if (obj.nStreams == 1) && nAntennas==1
                % in case of a single transmit antenna, no precoding is
                % needed
                if ~isequal( size(obj.precodingMatrix), [1, 1] )
                    obj.precodingMatrix = 1;
                    warning('For a single transmit antenna, the precoding matrix was set to 1!');
                end
                precodedDataSymbols = dataSymbols;
                return; 
            end
            if size(dataSymbols,2) ~= obj.nStreams
                error('Input data symbols size does not fit number of streams');
            end
            if (~strcmp(obj.method,'OLSM') && ~strcmp(obj.method,'TxD') && ~isequal( size(obj.precodingMatrix), [nAntennas, obj.nStreams] ))
                error('Precoding matrix does not fit number of antennas and number of streams!');
            end
            
            switch obj.method
                case 'CLSM'
                    precodedDataSymbols = dataSymbols * obj.codebook.W{obj.precodingMatrixIndex}.';
                case'OLSM'
                    %precodedDataSymbols = zeros(size(dataSymbols,1),obj.nStreams);
                    l = 1:size(dataSymbols,1);
                    if size(obj.codebook.W,1)==4
                        k = mod(floor((l-1)/obj.nStreams),4)+1;
                    else
                        k=ones(1,size(dataSymbols,1));
                    end
  
                    p = mod(l-1,obj.nStreams)+1;
                    for i=1:size(dataSymbols,1)
                        
                      currentPrecodingMatrix = obj.codebook.W{k(i),obj.nStreams}*(obj.codebook.D{obj.nStreams}).^(p(i))*obj.codebook.U{obj.nStreams};
                      precodedDataSymbols(i,:) = dataSymbols(i,:)*currentPrecodingMatrix.';
                         
                    end
                case 'TxD'
                    
                      % dim nStreams x Time
                     if ~ismember(nAntennas, [2, 4, 8])
                         error('TxD only supports 2, 4 or 8 transmit antennas');
                     end
                     
                     dataSymbols = reshape(dataSymbols,nAntennas,[]);   
                     X = 1/sqrt(2)*obj.codebook.Z{log2(nAntennas)}*[real(dataSymbols);imag(dataSymbols)];
                     c = length(X(1,:));
                     precodedDataSymbols = zeros(nAntennas,nAntennas*c);
                     for i = 1:nAntennas
                        precodedDataSymbols(:,i:nAntennas:nAntennas*c-(nAntennas-i))=X((i-1)*nAntennas+1:i*nAntennas,:);
                     end
                     precodedDataSymbols = precodedDataSymbols.';
                otherwise
                    precodedDataSymbols = dataSymbols * obj.precodingMatrix.';
            end
            % precoding

        end
        
         function setPrecodingMatrix( obj,precodingMatrix)
           
              obj.precodingMatrix = precodingMatrix;
               
         end
         
        function setPrecodingMatrixIndex( obj,PMI)
           
              obj.precodingMatrixIndex = PMI;
              obj.precodingMatrix = obj.codebook.W{PMI};
              obj.nStreams = rank(obj.precodingMatrix);
               
        end
        
        function precoder =getPrecoder(obj)
            if strcmp(obj.method,'OLSM')
                precoder.W = obj.codebook.W(:,obj.nStreams);
                precoder.U = obj.codebook.U{obj.nStreams};
                precoder.D = obj.codebook.D{obj.nStreams};
            else
                precoder = struct();
            end
            
        end
            
        
        
%         function setPrecodingMatrix( obj, type, varargin )
%             % this method updates the MIMO precoding matrix
%             switch type
%                 case 'index'
%                     % get LTE precoding matrix
%                 case 'custom'
%                     obj.precodingMatrix = varargin{1};
%                 otherwise
%                     error('Precoding method unknown!');
%             end
%         end
        
        function setNStreams( obj, nStreams )
            % set the number of spatial streams for transmission
            if isnumeric(nStreams)
                if nStreams > 0
                    obj.nStreams = nStreams;
                else
                    error('Number of streams must be positive!');
                end
            else
                error('Number of streams must be numeric!');
            end
        end
    end
    
end

