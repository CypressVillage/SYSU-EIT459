classdef Modulator < handle
%% The Modulator Class (Modulator)
%   Authors
%     - ...
%     - ...
%     - ...
%   (c) 2016 Institute of Telecommunications, TU Wien.
%   www.nt.tuwien.ac.at 

%% The Modulator Properties

% Parameters
    properties
        ModulationFormat                % Modulation Format: {'PAM', 'QAM'}
        ModulationOrder                 % Modulation Order: {4, 16, 64}
        Waveform                        % Waveform: {'OFDM', 'f-OFDM', 'FBMC', 'UFMC', 'GFDM'}
        ChannelEstimationMethod         % the channel estimation method: {'Approximate-Perfect', 'PilotAided'}
        noisePowerEstimation            % if the noise power is estimated {true,false}
        EqualizerType                   % The channel equalization method: {'One-Tap', 'FullBlock'}        
        ReceiverTypeMIMO                % The MIMO receiver method: {'ZF', 'MMSE', 'ML', 'Sphere'}        
        PilotSymbols                    % tranmitted pilot symbols
        symbolSpreading                 % data symbol spreading method
        channelAtPilots                 % the channel at pilot positions
        Channel                         % the estimated channel
    end
    
% Functionality
    properties
        WaveformObject                  % object for modulation and demodulation
        SignalConstellation             % object for bit to symbol mapping and de-mapping
        ChannelEstimator                % channel estimation object
        text
    end    
    
%% The Modulator Functions

% Setup
methods
    function obj = Modulator(varargin) 
        obj.ModulationFormat        = varargin{1};  % alphabet
        obj.ModulationOrder         = varargin{2};
        obj.Waveform                = varargin{3};
        obj.symbolSpreading         = varargin{4};
        obj.noisePowerEstimation    = varargin{5};
        obj.SignalConstellation     = Modulation.SignalConstellation(varargin{2}, varargin{1}, 0);
        
        switch obj.Waveform
            case 'OFDM'
                obj.WaveformObject = Modulation.OFDM(varargin{6}, varargin{7}, varargin{8}, varargin{9}, varargin{10}, varargin{11}, varargin{12}, varargin{13});
            case 'f-OFDM'
                obj.WaveformObject = Modulation.FOFDM(varargin{6}, varargin{7}, varargin{8}, varargin{9}, varargin{10}, varargin{11}, varargin{12}, varargin{13}, varargin{14}, varargin{15}, varargin{16});
            case 'WOLA'
                obj.WaveformObject = Modulation.WOLA(varargin{6}, varargin{7}, varargin{8}, varargin{9}, varargin{10}, varargin{11}, varargin{12}, varargin{13}, varargin{14}, varargin{15});
            case 'FBMC'
                obj.WaveformObject = Modulation.FBMC(varargin{6}, varargin{7}, varargin{8}, varargin{9}, varargin{10}, varargin{11}, varargin{12}, varargin{13}, varargin{14}, varargin{15});
            case 'UFMC'
                obj.WaveformObject = Modulation.UFMC(varargin{6}, varargin{7}, varargin{8}, varargin{9}, varargin{10}, varargin{11}, varargin{12}, varargin{13}, varargin{14}, varargin{15});
            case 'GFDM'
                obj.WaveformObject = Modulation.GFDM(varargin{6}, varargin{7}, varargin{8}, varargin{9}, varargin{10}, varargin{11}, varargin{12}, varargin{13});
            otherwise
                error('Waveform format unknown!');
        end
    end
    
    function setMCS( obj, modulationFormat, modulationOrder, MUSTIdx )
        % set the modulation format, order and constellation for adaptive
        % modulation and coding
        if strcmp(obj.Waveform, 'FBMC')
            obj.ModulationFormat    = 'PAM';
            obj.ModulationOrder     = sqrt(modulationOrder);
            obj.SignalConstellation = Modulation.SignalConstellation(sqrt(modulationOrder), 'PAM', 0);
        else
            obj.ModulationFormat    = modulationFormat;
            obj.ModulationOrder     = modulationOrder;
            obj.SignalConstellation = Modulation.SignalConstellation(modulationOrder, modulationFormat, MUSTIdx);
        end
    end
end

methods
    function nCodedBits = getNrCodedBits(obj, nPilotSymbols, schedule)
        % Get how many bits (coded) are required at the input of the modulator
        nCodedBits = (nnz(schedule)-nPilotSymbols)*log2(obj.ModulationOrder); 
    end
    
    function spreadSymbols = spreadingTransform(obj, dataSymbols, nScheduledSubcarriers)
        % this function performs a spreading transformation on the data
        % symbols
        
        % check input
        if ~isvector(dataSymbols)
            error('Input data symbols must be vector');
        end
        
        % check if spreading transform is turned on
        if strcmp(obj.symbolSpreading{1},'none')
            spreadSymbols = dataSymbols;
            return;
        end
        
        % reshape input vector
        dataSymbols = reshape(dataSymbols, nScheduledSubcarriers, []);
        
        % perform spreading
        switch obj.symbolSpreading{1}
            case{'DFT'}
                spreadSymbols = zeros(nScheduledSubcarriers, size(dataSymbols,2));
                for ii = 1:size(dataSymbols,2)
                    spreadSymbols(:,ii) = 1/sqrt(nScheduledSubcarriers) * fft(dataSymbols(:,ii));
                end
            otherwise
                error('Symbol spreading method unkown!');
        end
        
        % reshape to vector again
        spreadSymbols = spreadSymbols(:);
    end

    function despreadSymbols = inverseSpreadingTransform(obj, dataSymbols, nScheduledSubcarriers)
        % this function performs a spreading transformation on the data
        % symbols
        
        % check input
        if ~isvector(dataSymbols)
            error('Input data symbols must be vector');
        end
        
        % check if spreading transform is turned on
        if strcmp(obj.symbolSpreading{1},'none')
            despreadSymbols = dataSymbols;
            return;
        end
        
        % reshape input vector
        dataSymbols = reshape(dataSymbols, nScheduledSubcarriers, []);
        
        % perform spreading
        switch obj.symbolSpreading{1}
            case{'DFT'}
                despreadSymbols = zeros(nScheduledSubcarriers, size(dataSymbols,2));
                for ii = 1:size(dataSymbols,2)
                    despreadSymbols(:,ii) = sqrt(nScheduledSubcarriers) * ifft(dataSymbols(:,ii));
                end
            otherwise
                error('Symbol spreading method unkown!');
        end
        
        % reshape to vector again
        despreadSymbols = despreadSymbols(:);
    end
    
    function [LLRs,channelEstimationMSE] = demodulate(obj, channelOutput, schedule, noisePowerPerfect, MIMOmethod, precoder, precodingMatrix, channelObject, pathloss, txPower, powerScale)
        % This function calculates LLR values per bit from the received
        % signal. The obtained LLR values serve as an input for the soft 
        % channel decoder. It consists of demodulation, channel estimation,
        % equalzation and LLR calculation.
        
        % demodulate for each receive antenna individually
        nRx             = size(channelOutput,2);
        nTx             = size(precodingMatrix,1);
        nStreams        = size(precodingMatrix,2);
        pseudoSymbols   = zeros( obj.WaveformObject.Nr.Subcarriers, obj.WaveformObject.Nr.MCSymbols, nRx );
        for iRx = 1:nRx
            pseudoSymbols(:,:,iRx) = obj.WaveformObject.Demodulation(channelOutput(:,iRx));
        end
        
        % noise power estimation
        if obj.noisePowerEstimation
            interferenceNoisePower = zeros(nRx,obj.WaveformObject.Nr.Subcarriers,obj.WaveformObject.Nr.MCSymbols);
            for iRx = 1:nRx
                receivedSignalPadded = [channelOutput(:,iRx); zeros(obj.WaveformObject.Nr.MCSymbols-mod(size(channelOutput,1),obj.WaveformObject.Nr.MCSymbols),1)];
                % calculate received spectrum
                spectrum = circshift(1/sqrt(size(channelOutput,1)/obj.WaveformObject.Nr.MCSymbols) * fft(reshape(receivedSignalPadded,[],obj.WaveformObject.Nr.MCSymbols),obj.WaveformObject.Implementation.FFTSize),[-obj.WaveformObject.Implementation.IntermediateFrequency 0]);
                % get mean noise plus interference power from adjacent subcarriers
                interferenceNoisePower(iRx,:,:) = repmat(mean(abs([spectrum(end-(3:5),:); spectrum(obj.WaveformObject.Nr.Subcarriers+(3:5),:)]).^2,1),obj.WaveformObject.Nr.Subcarriers,1);
            end
            noisePower = interferenceNoisePower;
        else
            noisePower = noisePowerPerfect * ones(obj.WaveformObject.Nr.Subcarriers,obj.WaveformObject.Nr.MCSymbols);
        end
        
        % calculate perfect channel knowledge
        perfectChannel = 10^((txPower-30-pathloss)/20) * powerScale * channelObject.GetTransferFunction(obj.WaveformObject.GetTimeIndexMidPos,...
                                                                                                        obj.WaveformObject.Implementation.FFTSize,...
                                                                                                        obj.WaveformObject.Implementation.IntermediateFrequency+(1:obj.WaveformObject.Nr.Subcarriers));
        
        % equalization
        switch obj.EqualizerType
            case 'One-Tap'
                if strcmpi(obj.ChannelEstimationMethod,'Approximate-Perfect')
                    % get approximately perfect channel knowledge
                    obj.Channel = perfectChannel;
                    pilotMatrix = obj.ChannelEstimator.PilotMatrix;
                    dataIndex = ~sum(pilotMatrix,3);
                    channelEstimationMSE = NaN;
                    
                elseif strcmpi(obj.ChannelEstimationMethod,'PilotAided')
                     % LS channel estimation for each MIMO channel
                    pilotMatrix = obj.ChannelEstimator.PilotMatrix;
                    pilotsLS    = cell(nRx,nTx);
                    obj.Channel = zeros( size(pilotMatrix,1),  size(pilotMatrix,2), nRx, nTx );
                    for iRx = 1:nRx
                        pseudoSymbolsPerAnt = pseudoSymbols(:,:,iRx);
                        for iTx = 1:nTx
                            pilotsLS{iRx,iTx}       = zeros( obj.ChannelEstimator.NrPilotSymbols(iTx), 1 );
                            pilotsLS{iRx,iTx}       = pseudoSymbolsPerAnt(logical(pilotMatrix(:,:,iTx))) ./ obj.PilotSymbols{iTx};
                            obj.Channel(:,:,iRx,iTx)    = obj.ChannelEstimator.ChannelInterpolation(pilotsLS{iRx,iTx}, size(pilotMatrix,1), size(pilotMatrix,2), iTx);
                        end
                    end
                    dataIndex = ~sum(pilotMatrix,3);
                    
                    % calculate channel estimation error
                    channelEstimationMSE = 10^(-(txPower-30-pathloss)/10) * mean(abs(obj.Channel(:) - perfectChannel(:)).^2);
                    
                        % disp(obj.Channel)
                        % disp(size(obj.Channel))
                        x = 1:1:72;
                        y = 1:1:14;
                        [x,y] = meshgrid(x,y);
                        disp(size(x))
                        disp(size(y))
                        z = obj.Channel;
                        surf(x,y,10*log(abs(z')));
                        xlabel('子载波')
                        ylabel('OFDM符号')
                        zlabel('10*log_{10}|channel estimate|')
                        title('信道估计的频域响应')

                else
                    error('Channel estimation method not supported');
                end
                
                % get channel and demodulated data symbols (pseudo symbols)
                % at positions were data is scheduled
                
                scheduledChannel        = reshape( obj.Channel( repmat( dataIndex, 1, 1, channelObject.Nr.rxAntennas, channelObject.Nr.txAntennas ) ), [],  channelObject.Nr.rxAntennas, channelObject.Nr.txAntennas );
                scheduledPseudoymbos    = reshape( pseudoSymbols( repmat( dataIndex, 1,1, size(pseudoSymbols,3) ) ), [], size(pseudoSymbols,3) );
                obj.channelAtPilots=reshape( obj.Channel( repmat( ~dataIndex, 1, 1, channelObject.Nr.rxAntennas, channelObject.Nr.txAntennas ) ), [],  channelObject.Nr.rxAntennas, channelObject.Nr.txAntennas );
                % get some numbers
                nDataPositions          = nnz(dataIndex);
                nSymbols                = size(obj.Channel,2);
                nScheduledSubcarriers   = numel(scheduledChannel)/(nSymbols*nRx*nTx);
   
                % perform LLR calculation (detection)
                if (nTx==1) && (nRx==1)
                    % SISO detection => Zero forcing is optimal 
                    equalizedSymbols            = scheduledPseudoymbos./scheduledChannel;
                    symbolNoisePower            = obj.WaveformObject.GetSymbolNoisePower(squeeze(noisePower));
                    scheduledSymbolNoisePower   = reshape( symbolNoisePower( repmat( dataIndex, 1, 1, channelObject.Nr.rxAntennas, channelObject.Nr.txAntennas ) ), [],  channelObject.Nr.rxAntennas, channelObject.Nr.txAntennas );
                    scaledSymbolNoisePower      = scheduledSymbolNoisePower./abs(scheduledChannel).^2;
                    
                    % data symbol de-spreading
                    equalizedSymbolsDespread    = obj.inverseSpreadingTransform(equalizedSymbols, nScheduledSubcarriers);
                    
                    % LLR calculation
                    obj.text = equalizedSymbolsDespread(:);
                    LLRs = obj.SignalConstellation.LLR_AWGN( equalizedSymbolsDespread(:), scaledSymbolNoisePower(:) );
                else
                    % MIMO detection
                    switch obj.ReceiverTypeMIMO
                        case 'ZF'
                            LLRs = obj.SignalConstellation.LLR_MIMO_ZF(...
                                permute(scheduledPseudoymbos,[2 1]),...
                                permute(scheduledChannel,[2 3 1]),...
                                mean(mean(obj.WaveformObject.GetSymbolNoisePower(noisePower))),...
                                MIMOmethod, precoder, precodingMatrix,...
                                @obj.inverseSpreadingTransform,...
                                nScheduledSubcarriers);
                        case 'MMSE'
                            LLRs = obj.SignalConstellation.LLR_MIMO_MMSE(...
                                permute(scheduledPseudoymbos,[2 1]),...
                                permute(scheduledChannel,[2 3 1]),...
                                mean(mean(obj.WaveformObject.GetSymbolNoisePower(noisePower))),...
                                MIMOmethod, precoder, precodingMatrix,...
                                @obj.inverseSpreadingTransform,...
                                nScheduledSubcarriers);  
                        case 'Sphere'
                            if ~strcmp(obj.symbolSpreading,'none')
                                error('Currently shpere decoding is not supported in combination with data symbol spreading!');
                            end
                            LLRs = obj.SignalConstellation.LLR_MIMO_Sphere(...
                                permute(scheduledPseudoymbos,[2 1]),...
                                permute(scheduledChannel,[2 3 1]),...
                                mean(mean(obj.WaveformObject.GetSymbolNoisePower(noisePower))),...
                                MIMOmethod, precoder, precodingMatrix);
                        case 'ML'
                            if ~strcmp(obj.symbolSpreading,'none')
                                error('Currently maximum likelihood decoding is not supported in combination with data symbol spreading!');
                            end
                            LLRs = obj.SignalConstellation.LLR_MIMO_ML(...
                                permute(scheduledPseudoymbos,[2 1]),...
                                permute(scheduledChannel,[2 3 1]),...
                                repmat(eye(nRx),1,1,nDataPositions)*mean(mean(obj.WaveformObject.GetSymbolNoisePower(noisePower))),...
                                MIMOmethod, precoder, precodingMatrix);                         
                        otherwise
                              error('Equalizer type not supported');
                    end
                end
                          
            otherwise
                error('Equalizer type not supported');
        end  
    end
        
end
end 


