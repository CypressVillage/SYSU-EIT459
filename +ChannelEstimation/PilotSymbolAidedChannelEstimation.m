classdef PilotSymbolAidedChannelEstimation < handle
    % Ronald Nissel, rnissel@nt.tuwien.ac.at
    % (c) 2016 by Institute of Telecommunications, TU Wien
    % www.tc.tuwien.ac.at
    
    
    properties (SetAccess = private)
        NrPilotSymbols
        PilotPattern
        PilotSpacingFrequency
        PilotSpacingTime
        InterpolationMethod
        Implementation
        InterpolationProperties
        PilotMatrix
    end
    
    methods
        %% Class constructor, define default values.
        function obj = PilotSymbolAidedChannelEstimation(varargin)
            %% Initialize parameters, set default values
            obj.PilotPattern          = varargin{1};
            obj.InterpolationMethod   = varargin{3};
            schedule                  = varargin{4};
            nAntennas                 = varargin{5};
            
            % Generate pilot matrix according to the specified pilot pattern.
            % A zero corresponse to a data symbol, a one to a pilot symbol
            switch obj.PilotPattern
                case 'Rectangular'
                    nSubcarriers                = varargin{2}(1,1);
                    obj.PilotSpacingFrequency   = varargin{2}(1,2);
                    nMCSymbols                  = varargin{2}(2,1);
                    obj.PilotSpacingTime        = varargin{2}(2,2);
                    
                    % check number of antennas
                    if nAntennas > obj.PilotSpacingFrequency
                        error('Number of antennas must be smaller than frequency domain pilot spacing for the rectangular pilot pattern!');
                    end
                    
                    tmpPilotMatrix = zeros(nSubcarriers, nMCSymbols);
                    tmpPilotMatrix(round(mod(nSubcarriers-1,obj.PilotSpacingFrequency)/2)+1:obj.PilotSpacingFrequency:nSubcarriers,round(round(mod(nMCSymbols-1,obj.PilotSpacingTime)/2)+1:obj.PilotSpacingTime:nMCSymbols)) = true;
                    
                    obj.PilotMatrix   = zeros( nSubcarriers, nMCSymbols, nAntennas );
                    for iAnt = 1:nAntennas
                        obj.PilotMatrix(:,:,iAnt) = circshift(tmpPilotMatrix, iAnt-1, 1);
                    end
                case 'Diamond'
                    nSubcarriers = varargin{2}(1,1);
                    obj.PilotSpacingFrequency   = varargin{2}(1,2);
                    nMCSymbols                  = varargin{2}(2,1);
                    obj.PilotSpacingTime        = varargin{2}(2,2);
                    
                    % check number of antennas
                    if nAntennas > obj.PilotSpacingFrequency
                        error('Number of antennas must be smaller than frequency domain pilot spacing for the diamond pilot pattern!');
                    end
                                        
                    % There should be a much smarter way of doing this
                    tmpPilotMatrix = zeros(nSubcarriers, nMCSymbols);
                    FrequencyPositionShift = round((nSubcarriers-max([(1:2*obj.PilotSpacingFrequency:nSubcarriers),(1+1/2*obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:nSubcarriers),(1+obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:nSubcarriers),(1+3/2*obj.PilotSpacingFrequency:2*obj.PilotSpacingFrequency:nSubcarriers)]))/2)+1;
                    TimePositionShift = round((nMCSymbols-max([(1:2*obj.PilotSpacingTime:nMCSymbols),(1+obj.PilotSpacingTime):2*obj.PilotSpacingTime:nMCSymbols]))/2)+1;
                    tmpPilotMatrix(FrequencyPositionShift:2*obj.PilotSpacingFrequency:nSubcarriers,TimePositionShift:2*obj.PilotSpacingTime:nMCSymbols) = 1;
                    tmpPilotMatrix(FrequencyPositionShift+round(1/2*obj.PilotSpacingFrequency):2*obj.PilotSpacingFrequency:nSubcarriers,round(TimePositionShift+obj.PilotSpacingTime):2*obj.PilotSpacingTime:nMCSymbols) = 1;
                    tmpPilotMatrix(FrequencyPositionShift+round(obj.PilotSpacingFrequency):2*obj.PilotSpacingFrequency:nSubcarriers,TimePositionShift:2*obj.PilotSpacingTime:nMCSymbols) = 1;
                    tmpPilotMatrix(FrequencyPositionShift+round(3/2*obj.PilotSpacingFrequency):2*obj.PilotSpacingFrequency:nSubcarriers,round(TimePositionShift+obj.PilotSpacingTime):2*obj.PilotSpacingTime:nMCSymbols) = 1;
                    
                    obj.PilotMatrix   = zeros( nSubcarriers, nMCSymbols, nAntennas );
                    for iAnt = 1:nAntennas
                        obj.PilotMatrix(:,:,iAnt) = circshift(tmpPilotMatrix, iAnt-1, 1);
                    end
                case 'LTE Downlink'
                    nSubcarriers  = varargin{2}(1,1);
                    nMCSymbols    = varargin{2}(2,1);
                    
                    obj.PilotMatrix   = zeros( nSubcarriers, nMCSymbols, nAntennas );
                    pilotMatrixRb     = zeros( 12, nMCSymbols, nAntennas );
                    
                    % construct pilot pattern for antenna ports 1 and 2
                    unitCell = zeros(6,nMCSymbols/2);
                    unitCell(6,1)                 = 1;
                    unitCell(3,nMCSymbols/2-2)    = 1;
                    resourceBlock12 = repmat( unitCell, 2, 2);
                    
                    pilotMatrixRb(:,:,1)              = resourceBlock12;                          % antenna port 1
                    pilotMatrixRb(3:end,1:end-2,2)    = flipud( resourceBlock12(3:end,1:end-2) ); % antenna port 2
                    
                    % construct pilot pattern for antenna ports 3 and 4
                    resourceBlock3  = zeros(12,nMCSymbols);
                    resourceBlock3(6,2)               = 1;
                    resourceBlock3(3,nMCSymbols/2+2)  = 1;
                    resourceBlock3(12,2)              = 1;
                    resourceBlock3(9,nMCSymbols/2+2)  = 1;
                    
                    pilotMatrixRb(:,:,3)              = resourceBlock3;                           % antenna port 3
                    pilotMatrixRb(3:end,1:end-2,4)    = flipud( resourceBlock3(3:end,1:end-2) );  % antenna port 4
                    
                    % select antenna ports and repeat for each resource block
                    obj.PilotMatrix = repmat( pilotMatrixRb(:,:,1:nAntennas), nSubcarriers/12, 1, 1 );
                    
                case 'Custom'
                    obj.PilotSpacingFrequency = nan;
                    obj.PilotSpacingTime      = nan;
                    obj.PilotMatrix           = varargin{2};
                otherwise
                    error('Pilot pattern is not supported!');
            end
%             obj.PilotMatrix = logical( obj.PilotMatrix & repmat( schedule, 1, 1, nAntennas ) );
            obj.NrPilotSymbols = squeeze(sum(sum(obj.PilotMatrix,1),2));
            
%             % pre initialize interpolator
%             % this should be done per antenna for MIMO
%             switch obj.InterpolationMethod
%                 case {'linear','nearest','natural'}
%                     [x_pilot_pos,y_pilot_pos] = find(obj.PilotMatrix);
%                     obj.InterpolationProperties = scatteredInterpolant(x_pilot_pos,y_pilot_pos,zeros(mean(obj.NrPilotSymbols),1),obj.InterpolationMethod);
%                 case 'MovingBlockAverage'
%                     PilotIndices = find(obj.PilotMatrix);
%                     PilotMatrixPosNumbered = zeros(size(obj.PilotMatrix));
%                     PilotMatrixPosNumbered(PilotIndices)=1:numel(PilotIndices);
%                     
%                     BlockLengthFrequency = varargin{4}(1);
%                     BlockLengthTime = varargin{4}(2);
%                     
%                     maxF = size(obj.PilotMatrix,1);
%                     maxT = size(obj.PilotMatrix,2);
%                     
%                     IndexF = -BlockLengthFrequency:BlockLengthFrequency;
%                     IndexT = -BlockLengthTime:BlockLengthTime;
%                     
%                     obj.InterpolationProperties.InterpolationMatrix = zeros(numel(obj.PilotMatrix),obj.NrPilotSymbols);
%                     for i_pos = 1:numel(obj.PilotMatrix)
%                         Impulse = zeros(size(obj.PilotMatrix));
%                         Impulse(i_pos)=1;
%                         [posF,posT]=find(Impulse);
%                         
%                         IndexPosT = posT+IndexT;
%                         IndexPosT(IndexPosT<1)=[];
%                         IndexPosT(IndexPosT>maxT)=[];
%                         IndexPosF = posF+IndexF;
%                         IndexPosF(IndexPosF<1)=[];
%                         IndexPosF(IndexPosF>maxF)=[];
%                         Impulse(IndexPosF,IndexPosT)=1;
%                         
%                         InterpolationMatrixPosPilots = PilotMatrixPosNumbered(logical(Impulse) & logical(obj.PilotMatrix));
%                         
%                         obj.InterpolationProperties.InterpolationMatrix(i_pos,InterpolationMatrixPosPilots) = 1/numel(InterpolationMatrixPosPilots);
%                     end
%                 case 'MMSE'
%                     error('Needs to be done');
%             end
        end
        
        function InterpolatedChannel = ChannelInterpolation(obj, LSChannelEstimatesAtPilotPosition, BlockLengthFrequency, BlockLengthTime, iAntenna)
            % interpolation method
            switch obj.InterpolationMethod
                case {'linear','nearest','natural'}
                    [x_pilot_pos,y_pilot_pos] = find(obj.PilotMatrix(:,:,iAntenna));
                    obj.InterpolationProperties = scatteredInterpolant(x_pilot_pos,y_pilot_pos,zeros(obj.NrPilotSymbols(iAntenna),1),obj.InterpolationMethod);
                case 'MovingBlockAverage'
                    PilotIndices = find(obj.PilotMatrix(:,:,iAntenna));
                    PilotMatrixPosNumbered = zeros(size(obj.PilotMatrix(:,:,iAntenna)));
                    PilotMatrixPosNumbered(PilotIndices)=1:numel(PilotIndices);
                    
                    maxF = size(obj.PilotMatrix,1);
                    maxT = size(obj.PilotMatrix,2);
                    
                    IndexF = -BlockLengthFrequency:BlockLengthFrequency;
                    IndexT = -BlockLengthTime:BlockLengthTime;
                    
                    obj.InterpolationProperties.InterpolationMatrix = zeros(numel(obj.PilotMatrix),obj.NrPilotSymbols(iAntenna));
                    for i_pos = 1:numel(obj.PilotMatrix(:,:,iAntenna))
                        Impulse = zeros(size(obj.PilotMatrix(:,:,iAntenna)));
                        Impulse(i_pos)=1;
                        [posF,posT]=find(Impulse);
                        
                        IndexPosT = posT+IndexT;
                        IndexPosT(IndexPosT<1)=[];
                        IndexPosT(IndexPosT>maxT)=[];
                        IndexPosF = posF+IndexF;
                        IndexPosF(IndexPosF<1)=[];
                        IndexPosF(IndexPosF>maxF)=[];
                        Impulse(IndexPosF,IndexPosT)=1;
                        
                        InterpolationMatrixPosPilots = PilotMatrixPosNumbered(logical(Impulse) & logical(obj.PilotMatrix(:,:,iAntenna)));
                        
                        obj.InterpolationProperties.InterpolationMatrix(i_pos,InterpolationMatrixPosPilots) = 1/numel(InterpolationMatrixPosPilots);
                    end
                otherwise
                    error('Channel interpolation method not supported!');
            end
            
            switch obj.InterpolationMethod
                case {'linear','nearest','natural'}
                    obj.InterpolationProperties.Values = LSChannelEstimatesAtPilotPosition;
                    [yq,xq] = meshgrid(1:size(obj.PilotMatrix,2),1:size(obj.PilotMatrix,1));
                    InterpolatedChannel = obj.InterpolationProperties(xq,yq);
                case 'FullAverage'
                    InterpolatedChannel = ones(size(obj.PilotMatrix))*mean(LSChannelEstimatesAtPilotPosition);
                case 'MovingBlockAverage'
                    InterpolatedChannel = obj.InterpolationProperties.InterpolationMatrix*LSChannelEstimatesAtPilotPosition;
                otherwise
                    error('Interpolation method not implemented');
            end
        end
        
        
        
        function AuxiliaryMatrix = GetAuxiliaryMatrix(varargin)
            obj = varargin{1};
            NrAxuiliarySymbols = varargin{2};
            
            AuxiliaryMatrix = obj.PilotMatrix;
            [index_l,index_k]=find(obj.PilotMatrix);
            if (min(index_l)<2) || (max(index_l)>=size(obj.PilotMatrix,1))
                warning('Pilots should not be close to the border! There might be a problem!');
            elseif (min(index_k)<2) || (max(index_k)>=size(obj.PilotMatrix,2))
                warning('Pilots should not be close to the border! There might be a problem!');
            end
            
            for i_lk = 1:size(index_l,1)
                switch NrAxuiliarySymbols
                    case 1
                        AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)+1) = -1;
                    case 2
                        AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)+1) = -1;
                        AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)-1) = -1;
                    case 3
                        AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)+1) = -1;
                        AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)-1) = -1;
                        AuxiliaryMatrix(index_l(i_lk)+1,index_k(i_lk)) = -1;
                    case 4
                        AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)+1) = -1;
                        AuxiliaryMatrix(index_l(i_lk),index_k(i_lk)-1) = -1;
                        AuxiliaryMatrix(index_l(i_lk)+1,index_k(i_lk)) = -1;
                        AuxiliaryMatrix(index_l(i_lk)-1,index_k(i_lk)) = -1;
                    otherwise
                        error('Only 1,2,3,4 auxiliary symbols per pilot are supported');
                end
            end
        end
        
        function InterpolationMatrix = GetInterpolationMatrix(varargin)
            obj = varargin{1};
            
            [x_pilot_pos,y_pilot_pos] = find(obj.PilotMatrix);
            InterpolationMatrix = zeros(numel(obj.PilotMatrix),numel(x_pilot_pos));
            for i_pos =1:length(x_pilot_pos)
                TestDirac = zeros(size(x_pilot_pos));
                TestDirac(i_pos)=1;
                ImpulseResponse = obj.ChannelInterpolation(TestDirac);
                
                InterpolationMatrix(:,i_pos)=ImpulseResponse(:);
            end
            
        end
        
        function PlotPilotPattern(varargin)
            if numel(varargin)==2
                PilotMatrixTemp = varargin{2};
            else
                PilotMatrixTemp = varargin{1}.PilotMatrix;
            end
            
            PilotMatrixRGB(:,:,1) = PilotMatrixTemp==0;
            PilotMatrixRGB(:,:,2) = not(PilotMatrixTemp==1);
            PilotMatrixRGB(:,:,3) = not(PilotMatrixTemp==-1);
            
            imagesc(PilotMatrixRGB);
            hold on;
            for i_row = 1:size(PilotMatrixRGB,1)+1
                plot([.5,size(PilotMatrixRGB,2)+0.5],[i_row-.5,i_row-.5],'k-');
            end
            for i_column = 1:size(PilotMatrixRGB,2)+1
                plot([i_column-.5,i_column-.5],[.5,size(PilotMatrixRGB,1)+0.5],'k-');
            end
            xlabel('Time index');
            ylabel('Frequency index');
            title('Pilot pattern');
        end
        
    end
    
    
end