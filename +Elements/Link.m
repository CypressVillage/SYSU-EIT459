classdef Link < handle
    %% The Link Class (Link)
    %   Authors
    %     - ...
    %     - ...
    %     - ...
    %   (c) 2016 Institute of Telecommunications, TU Wien.
    %   www.nt.tuwien.ac.at
    
    %% The Link Properties
    % The properites are divided into five categories:
    
    % Accessibility
    properties
        ID                      % Link ID.
        Type                    % Link Type: {'Primary', 'Interference'}.
        direction               % 'downlink','uplink','D2D'
        Transmitter             % The transmitting end (node ID).
        Receiver                % The receiving end (node ID).
        MUSTFarUE               % ID of the MUST FarUE.
        MUSTNearUE              % ID of the MUST NearUE.
    end
    
    % Parameters
    properties
        SNR                     % The SNR at the receiving end
        pathloss                % artificial channel pathloss
        transmissionMode        % 'CLSM', 'OLSM' or TxD
        HARQType                % The HARQ type: {'CC', 'IR'}.
        HARQRV = 0;             % The HARQ Redudency Version: {0, 1, 2, ..., max(HARQRV)}
        Velocity = 0;           % The relative velocity between the transmitting and receiving ends
        ChannelEstimationMethod % The channel estimation method: {'Perfect', 'PilotAided'}
        NrInputBits             % Number of the original data bits
        schedule                % schedule for transmission in current frame
        scheduledSubcarriers    % number of total scheduled subcarriers
        firstSubcarrier         % position of the first scheduled subcarrier, ie., the frequency offset in subcarriers
        isScheduled = 0         % Indicates whether the links is Scheduled or not: isScheduled = {0, 1};
        Attenuation = 0         % Scales down the link output power by a factor (in dB). ex: Attenuation = 5;
        transmitPower           % link's transmit power
        powerScale              % power scaling to get a total power of one for all schedules (for downlink)
        nTxAntennas             % number of transmit antennas
        nRxAntennas             % number of receive antennas
        IsMUSTFarUE = 0;
        IsMUSTNearUE = 0;
        layerMappingTable
        Nonlinearity            % Apply nonlinear model
        smoothnessFactor        % Smoothness factor of the model
        amplifierOBO            % Amplifier Output Back-off

    end
    
    % Functionality
    properties
        ChannelCoder            % The channel coding object.
        Modulator               % The modulation object.
        MIMO                    % The MIMO object.
        Channel                 % The channel between the transmitting and receiving end.
        Feedback                % The Feedback object
    end
    
    % Generated Signals
    properties
        InputBits
        EncodedBits
        TransmitSymbols
        TransmitSignal
        ReceiveSignal
        LLRs
        DecodedBits
        channelEstimationMSE
    end
    
    % Results of current transmission
    properties
        BER                % BER
        FrameError         % FE: {0, 1}
        Throughput         % Throughput
    end
    %% The Link Functions
    % The Link provide functions for both the transmit and receive sides
    
    % Setup
    methods

        function obj = Link(ID, type, direction, transmitter, receiver, harqType, velocity, attenuation,  nTXAntennas, nRXAntennas, feedback, modulation, transmissionMode)
            obj.ID              = ID;
            obj.Type            = type;
            obj.direction       = direction;
            obj.Transmitter     = transmitter;
            obj.Receiver        = receiver;
            obj.HARQType        = harqType;
            obj.Velocity        = velocity;
            obj.Attenuation 	= attenuation;
            obj.nTxAntennas     = nTXAntennas;
            obj.nRxAntennas     = nRXAntennas;
            obj.transmissionMode = transmissionMode;
            if  (~strcmp(type,'Interference') && feedback.enable)
                obj.Feedback        = Feedback.Feedback(feedback, nTXAntennas, nRXAntennas,  modulation.precodingMatrix{ID}); 
            end
        end
        
        function updateLink(obj, simParams, Links, frame)
            % In this method, the current link is updated: scheduling and
            % link adaptation, feedback calculation
            
            %% Setting of the Channel Estimator and Equalizer
            obj.Modulator.EqualizerType = simParams.simulation.equalizerType;
            obj.Modulator.ReceiverTypeMIMO = simParams.simulation.receiverTypeMIMO;            
            obj.ChannelEstimationMethod = simParams.simulation.channelEstimationMethod;
            switch obj.ChannelEstimationMethod
                case 'Approximate-Perfect'
                    % artificially add pilots to get the same throughput as
                    % with channel estimation
                    ChannelEstimator = ChannelEstimation.PilotSymbolAidedChannelEstimation( ...
                        simParams.simulation.pilotPattern,...           % Pilot pattern
                        [...                                            % Matrix that represents the pilot pattern parameters
                        obj.Modulator.WaveformObject.Nr.Subcarriers,... % Number of subcarriers
                        6; ...                                          % Pilot spacing in the frequency domain
                        obj.Modulator.WaveformObject.Nr.MCSymbols,...   % Number of Symbols
                        3.5 ...                                         % Pilot spacing in the time domain
                        ],...
                        'linear',...                                    % Interpolation(Extrapolation) method
                        obj.schedule, ...
                        obj.nTxAntennas ...
                        );
                     
                    obj.Modulator.ChannelEstimationMethod = 'Approximate-Perfect';
                    obj.Modulator.ChannelEstimator = ChannelEstimator;
                    nPilotSymbols = ChannelEstimator.NrPilotSymbols;               
                case 'Perfect'
                    obj.Modulator.ChannelEstimationMethod = 'Perfect';
                    nPilotSymbols = 0;
                case 'PilotAided'
                    ChannelEstimator = ChannelEstimation.PilotSymbolAidedChannelEstimation( ...
                        simParams.simulation.pilotPattern,...           % Pilot pattern
                        [...                                            % Matrix that represents the pilot pattern parameters
                        obj.Modulator.WaveformObject.Nr.Subcarriers,... % Number of subcarriers
                        1; ...                                          % Pilot spacing in the frequency domain
                        obj.Modulator.WaveformObject.Nr.MCSymbols,...   % Number of Symbols
                        6 ...                                         % Pilot spacing in the time domain
                        ],...
                        'linear',...                                    % Interpolation(Extrapolation) method
                        obj.schedule, ...
                        obj.nTxAntennas ...
                        );
                    obj.Modulator.ChannelEstimationMethod = 'PilotAided';
                    obj.Modulator.ChannelEstimator = ChannelEstimator;
                    nPilotSymbols = ChannelEstimator.NrPilotSymbols;

                    ChannelEstimator.PlotPilotPattern();
                otherwise
                    error('Channel estimation method unknown.');
            end
            
            %% Feedback
            obj.Channel.NewRealization(frame);
            switch obj.transmissionMode
                case 'CLSM'
                    % feedback is set based on the transmission mode. 
                    if (~isempty(obj.powerScale))
                        % if the feedback delay is 0, the channel for the next
                        % transmission is used for the feedback calculation else
                        % the current channel is used
                        if obj.Feedback.delay==0
                            channel = 10^((obj.transmitPower-30-obj.pathloss)/20) * obj.powerScale * obj.Channel.GetTransferFunction(obj.Modulator.WaveformObject.GetTimeIndexMidPos,...
                                obj.Modulator.WaveformObject.Implementation.FFTSize,...
                                obj.Modulator.WaveformObject.Implementation.IntermediateFrequency+(1:obj.Modulator.WaveformObject.Nr.Subcarriers));  
                        else
                            channel=obj.Modulator.Channel;
                        end
                        if obj.IsMUSTNearUE
                                channel = channel*obj.Modulator.SignalConstellation.MUSTNearUEScale;
                        end
                        pilotMatrix = obj.Modulator.ChannelEstimator.PilotMatrix; 
                        noisePower = obj.Modulator.WaveformObject.GetSymbolNoisePower(simParams.phy.noisePower);
                        % instead of noisePower simParams.phy.noisePower
                        obj.Feedback.updateFeedback(channel,pilotMatrix,noisePower(1),[simParams.modulation.mcsValues.modulationOrder]);
                    end
                        mcs = simParams.modulation.mcsValues(obj.Feedback.Cqi.cqiArray{1}(1));
                        obj.MIMO.setNStreams(obj.Feedback.Ri.riArray(1));
                        obj.MIMO.setPrecodingMatrixIndex( obj.Feedback.Pmi.pmiArray(1) );
                                        
                case 'OLSM'
                    if (~isempty(obj.powerScale))

                        if obj.Feedback.delay==0
                            channel = 10^((obj.transmitPower-30-obj.pathloss)/20) * obj.powerScale * obj.Channel.GetTransferFunction(obj.Modulator.WaveformObject.GetTimeIndexMidPos,...
                                obj.Modulator.WaveformObject.Implementation.FFTSize,...
                                obj.Modulator.WaveformObject.Implementation.IntermediateFrequency+(1:obj.Modulator.WaveformObject.Nr.Subcarriers));  
                        else
                            channel=obj.Modulator.Channel;
                        end
                        if obj.IsMUSTNearUE
                                channel = channel*obj.Modulator.SignalConstellation.MUSTNearUEScale;
                        end

                        pilotMatrix = obj.Modulator.ChannelEstimator.PilotMatrix; 
                        noisePower = obj.Modulator.WaveformObject.GetSymbolNoisePower(simParams.phy.noisePower);
 
                        obj.Feedback.updateFeedback(channel,pilotMatrix,noisePower(1),[simParams.modulation.mcsValues.modulationOrder]);
                    end
                    
                        mcs = simParams.modulation.mcsValues(obj.Feedback.Cqi.cqiArray{1}(1)); 
                        obj.MIMO.setNStreams(obj.Feedback.Ri.riArray(1));
                        
                    
                case 'TxD'
                    
                    mcs = simParams.modulation.mcsValues(simParams.modulation.mcs(obj.ID) );
                    obj.MIMO.setNStreams(simParams.modulation.nStreams(obj.ID) );
                    obj.MIMO.setPrecodingMatrix(simParams.modulation.precodingMatrix{obj.ID} );
                
                otherwise
                    % feedback up date for the custom transmission mode
                        if (~isempty(obj.Feedback) && ~isempty(obj.powerScale))
                            if obj.Feedback.delay==0
                                channel = 10^((obj.transmitPower-30-obj.pathloss)/20) * obj.powerScale * obj.Channel.GetTransferFunction(obj.Modulator.WaveformObject.GetTimeIndexMidPos,...
                                    obj.Modulator.WaveformObject.Implementation.FFTSize,...
                                    obj.Modulator.WaveformObject.Implementation.IntermediateFrequency+(1:obj.Modulator.WaveformObject.Nr.Subcarriers));  
                            else
                                channel=obj.Modulator.Channel;
                            end

                            if obj.IsMUSTNearUE
                                channel = channel*obj.Modulator.SignalConstellation.MUSTNearUEScale;
                            end

                            pilotMatrix = obj.Modulator.ChannelEstimator.PilotMatrix; 

                            obj.Feedback.updateFeedback(channel,pilotMatrix,simParams.phy.noisePower,[simParams.modulation.mcsValues.modulationOrder]);

                            if obj.Feedback.Cqi.enable
                                mcs = simParams.modulation.mcsValues(obj.Feedback.Cqi.cqiArray{1}(1));
                            else
                                mcs = simParams.modulation.mcsValues(simParams.modulation.mcs(obj.ID));
                            end

                            if obj.Feedback.Pmi.enable
                                obj.MIMO.setNStreams(obj.Feedback.Ri.riArray(1));
                                obj.MIMO.setPrecodingMatrix( obj.Feedback.getPrecodingMatrix(obj.Feedback.Pmi.pmiArray(1)) );
                            else
                                obj.MIMO.setNStreams(simParams.modulation.nStreams(obj.ID));
                                obj.MIMO.setPrecodingMatrix(simParams.modulation.precodingMatrix{obj.ID});
                            end
                        else
                            % feedback deactivated 
                            mcs = simParams.modulation.mcsValues(simParams.modulation.mcs(obj.ID));
                            obj.MIMO.setNStreams(simParams.modulation.nStreams(obj.ID));
                            obj.MIMO.setPrecodingMatrix(simParams.modulation.precodingMatrix{obj.ID});
                        end
            end                

            %% MUST Operation
            if obj.IsMUSTFarUE
                % If UE is MUSTFarUE then the mapping is limited to QPSK (CQI 6) 
                if mcs.cqi > 6
                    mcs = simParams.modulation.mcsValues(6);
                end
                % Forces FarUE to use the same MIMO setting as NearUE.
                obj.MIMO = Links{obj.Transmitter, obj.MUSTNearUE}.MIMO;
            end
        
            % If MUST is enabled, choose the power ratio (affects NearUE only).
            if obj.IsMUSTNearUE
                % Basic choice for now, choose the middle ratio (non-adaptive).
                MUSTIdx = 2;
            else
                MUSTIdx = 0;
            end
            %% Link Adaptation
            obj.Modulator.setMCS('QAM', mcs.modulationOrder, MUSTIdx);
            nCodedBits = obj.Modulator.getNrCodedBits(sum(nPilotSymbols), obj.schedule);
            
            % Limited Soft Buffer: 1 for 100% buffer size (no limitation).
            SoftBufferRatio = 1;
            
            switch obj.direction
                case 'Uplink'
                     layerMapping = obj.layerMappingTable.Uplink{obj.MIMO.nStreams};
                case 'Downlink'
                    layerMapping = obj.layerMappingTable.Downlink{obj.MIMO.nStreams};
                otherwise
                    error('Layer mapping table unknown');
            end
            nCodewords = length(layerMapping);
            
            for iCodeword = 1:nCodewords 
                obj.NrInputBits(iCodeword) = obj.ChannelCoder{iCodeword}.update('Output', nCodedBits*layerMapping(iCodeword), mcs.codingRateTimes1024/1024, log2(mcs.modulationOrder), SoftBufferRatio);
            end
                         
        end
    end
    
    
    % Transmit and channel side
    methods
        function generateReceiveSignal(obj)    
            switch obj.Type
                case 'Primary'
                    obj.ReceiveSignal =                             obj.Channel.Convolution( (10^(-obj.pathloss/20)) * obj.TransmitSignal );
                case 'Interference'
%                     obj.ReceiveSignal = 10^(-obj.Attenuation/20) *  obj.Channel.Convolution( (10^(-obj.pathloss/20)) * obj.TransmitSignal );
                    obj.ReceiveSignal = 10^(-obj.Attenuation/20) *  obj.Channel.Convolution(obj.TransmitSignal ); % Needs further work!
            end
            
        end
        
        function applyNonlinearity(obj, transmitSignal)
            if obj.Nonlinearity
                transmitSignal = reshape(transmitSignal,[],1);                                      % Parallel handling of all nTx streams
                limitingAmplitude = 10.^((obj.amplifierOBO)/20);                                    % Amplifier saturation amplitude in linear units
                inputAmplitude = sqrt(real(transmitSignal).^2 + imag(transmitSignal).^2);           % Calculate each symbol's amplitude
                outputAmplitude = inputAmplitude ./ (1 + (inputAmplitude/limitingAmplitude).^(2*obj.smoothnessFactor)).^(1/(2*obj.smoothnessFactor));    %apply nonlinear model on the signal
                scalingFactor = outputAmplitude ./ inputAmplitude;                                  % calculate scaling factor and use it on the signal vector
                scalingFactor(isnan(scalingFactor)) = 0;                                            % ommit NAN's in case of division by 0 in previous line
                obj.TransmitSignal = reshape(transmitSignal .* scalingFactor,[],obj.nTxAntennas);   % distort signal vector
                
                meanPower = norm(obj.TransmitSignal)^2 / length(obj.TransmitSignal);                % calculate mean signal power
                obj.TransmitSignal = obj.TransmitSignal/sqrt(meanPower);                            % normalize output signal power to 1 (this seems like normalizing to the mean of all antennas)
                
            else                                                                                    % return output signal without alternation
                obj.TransmitSignal = transmitSignal;
            end
        end
    end
    
    % Receive-side
    methods
        function demodulate(obj, receivedSignal, noisePower)

            [obj.LLRs, obj.channelEstimationMSE] = obj.Modulator.demodulate(receivedSignal,...
                                                                            obj.schedule,...
                                                                            noisePower,...
                                                                            obj.MIMO.method,...
                                                                            obj.MIMO.getPrecoder(),...
                                                                            obj.MIMO.precodingMatrix,...
                                                                            obj.Channel,...
                                                                            obj.pathloss,...
                                                                            obj.transmitPower,...
                                                                            obj.powerScale); % The last three arguments are only necessary for perfect channel knowledge
                                                                        
            switch obj.direction
                case 'Uplink'
                     layerMapping = obj.layerMappingTable.Uplink{obj.MIMO.nStreams};
                case 'Downlink'
                    layerMapping = obj.layerMappingTable.Downlink{obj.MIMO.nStreams};
                otherwise
                    error('Layer mapping table unknown');
            end
                
            nCodewords = length(layerMapping);
            LLRsPerCodeword = cell(nCodewords,1);
            
            endLayer = 0;
            for iCodeword=1:nCodewords
                startLayer = endLayer+1;
                endLayer = startLayer+layerMapping(iCodeword)-1;
                nLayers = layerMapping(iCodeword);
                LLRsPerLayer = obj.LLRs(:,startLayer:endLayer);
                nTransmitSymbolsPerLayer = length(obj.TransmitSymbols{iCodeword});
                LLRsPerCodewordTemp = zeros(log2(obj.Modulator.SignalConstellation.ModulationOrder),nTransmitSymbolsPerLayer);
                for iLayer =1:nLayers
                    LLRsPerCodewordTemp(:,iLayer:nLayers:end) = reshape(LLRsPerLayer(:,iLayer),log2(obj.Modulator.SignalConstellation.ModulationOrder),[]);
                end
                LLRsPerCodeword{iCodeword} = reshape(LLRsPerCodewordTemp,1,[]);
                
                % If MUST is enabled (perform ML NOMA detection)
                if obj.IsMUSTNearUE
                    compositeLLRs = reshape(LLRsPerCodeword{iCodeword}', log2(obj.Modulator.ModulationOrder) + 2, [])';

                    % Remove FarUE bits
                    compositeLLRs = compositeLLRs(:, 1:log2(obj.Modulator.ModulationOrder));
                    
                    LLRsPerCodeword{iCodeword} = compositeLLRs(:).';
                end 
            end
            obj.LLRs = LLRsPerCodeword;       
        end
        
        function decode(obj)
            % decoding for each spatial stream individually
            nCodewords = length(obj.LLRs);
            obj.DecodedBits = cell(nCodewords,1);
            obj.FrameError  = zeros(nCodewords,1);
            
            for iCodeword = 1:nCodewords
                % decoding
                obj.DecodedBits{iCodeword} = obj.ChannelCoder{iCodeword}.decode(-obj.LLRs{iCodeword});
                % check if codeword was in error
                obj.FrameError(iCodeword) = obj.ChannelCoder{iCodeword}.CRCDetectionResult;
            end
        end
    end
   
    
    % Results
    methods
        function calculateSNR( obj, BOLTZMANN, temperature)
            % this function calculates the SNR
            if strcmp(obj.Type, 'Interference')
                error('There is no SNR calculation for an interference link.');
            end
            
            if strcmp(obj.direction, 'Downlink')
                % for downlink the total transmit power is spread over all
                % scheduled subcarriers (for all users combined)
                nSubcarriers = mean(obj.scheduledSubcarriers);
            elseif strcmp(obj.direction, 'Uplink') || strcmp(obj.direction, 'D2D')
                % for uplink the total transmit power is spread over all
                % scheduled subcarriers of the user
                nSubcarriers = mean(sum(double(obj.schedule)));
            else
                error('Link direction unkown');
            end

            noisePowerDensity   = BOLTZMANN * temperature;
            signalPowerDensity  = 10^( (obj.transmitPower - 30 - obj.pathloss)/10 ) / ( nSubcarriers * obj.Modulator.WaveformObject.PHY.SubcarrierSpacing );
            obj.SNR             = 10*log10( signalPowerDensity / noisePowerDensity );
        end
        
        function results = getResults(obj, saveData)
            % Return results
            results.SNR                 = obj.SNR;
            results.nStreams            = obj.MIMO.nStreams;
            switch obj.direction
                    case 'Uplink'
                         layerMapping = obj.layerMappingTable.Uplink{obj.MIMO.nStreams};
                    case 'Downlink'
                         layerMapping = obj.layerMappingTable.Downlink{obj.MIMO.nStreams};
                    otherwise
                         error('Layer mapping table unknown');
            end 
            nCodewords = length(layerMapping);
            results.nDataBits = obj.NrInputBits;
            results.frameError = obj.FrameError;
            results.channelMSE = obj.channelEstimationMSE;
            encodedBitsReceived = cell(nCodewords,1);
            for iCodeword = 1:nCodewords
                results.nCodeBits(iCodeword)= length( obj.EncodedBits{iCodeword} );
                encodedBitsReceived{iCodeword} = sign(obj.LLRs{iCodeword})==1;
                results.nBitError(iCodeword) = sum(obj.EncodedBits{iCodeword}~=encodedBitsReceived{iCodeword}.');
                results.nCodedBitError(iCodeword) = sum(obj.InputBits{iCodeword}~=obj.DecodedBits{iCodeword});
            end
            
            if saveData
                results.transmittedSymbols  = obj.TransmitSymbols;
                results.dataBits            = obj.InputBits;
                results.encodedBits         = obj.EncodedBits;
                results.encodedBitsReceived = encodedBitsReceived;
                results.decodedBits         = obj.DecodedBits;  
            end
            
            %PAPR calculation section
            oversamplingFactor = 4;                                          % minimal oversampling factor to ensure correct PAPR calculation
            timeVector = (0:obj.Modulator.WaveformObject.PHY.dt:(obj.Modulator.WaveformObject.PHY.dt*(obj.Modulator.WaveformObject.Nr.SamplesTotal  -1))).';
            timeVectorDesired = ((0:ceil(timeVector(end)/(obj.Modulator.WaveformObject.PHY.dt / oversamplingFactor)))*(obj.Modulator.WaveformObject.PHY.dt / oversamplingFactor)).';
            transmitSignalInterpolated = interp1(timeVector,obj.TransmitSignal,timeVectorDesired);
            results.PAPR = 10*log10(max(max(abs(transmitSignalInterpolated).^2)));
        end
    end
    
end

