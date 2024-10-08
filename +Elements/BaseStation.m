classdef BaseStation < handle & Elements.Node
%% The BaseStation Class (BaseStation)
%   Authors
%     - ...
%     - ...
%   (c) 2016 Institute of Telecommunications, TU Wien.
%   www.nt.tuwien.ac.at 
%% Basetation Properties
    properties
        txPower
        nAntennas
    end
%% The BaseStation Functions
%% Initialization
    methods
        function obj = BaseStation(name, id, nAntennas, transmitPower, plotResults)
            obj = obj@Elements.Node(name, id, 'BaseStation', plotResults);
            obj.nAntennas   = nAntennas;
            obj.txPower     = transmitPower;
        end
    end
%% Transmit Functions
    methods
        function generateTransmitSignal(obj, Links)
                        
            totalSubcarriers            = 0;
            totalScheduledSubcarriers   = 0;
            iScheduledUE                = 0;
            
            % Loop through all schduelled users and generate the total transmit signal
            for iUE = 1:length(obj.ReceiveUE)
                currentLink = Links{obj.ID, obj.ReceiveUE(iUE)};
                if any(currentLink.schedule(:)) && strcmp(currentLink.Type, 'Primary') && ~currentLink.IsMUSTFarUE
                    currentLink.isScheduled = 1;
                    iScheduledUE = iScheduledUE + 1;
                          
                    % get total an scheduled subcarriers
                    totalSubcarriers = size(currentLink.schedule,1);
                    totalScheduledSubcarriers = totalScheduledSubcarriers + currentLink.scheduledSubcarriers;
                    
                    nStreams = currentLink.MIMO.nStreams;
                    
                    % TODO different sized codewords
                    nCodewords = length(currentLink.layerMappingTable.Downlink{nStreams});
                    
                    % Channel Coding for this user
                    currentLink.InputBits   = cell(nCodewords,1);
                    currentLink.EncodedBits = cell(nCodewords,1);

                    for iCodeword = 1:nCodewords
                        currentLink.InputBits{iCodeword} =  randi([0 1], currentLink.NrInputBits(iCodeword), 1);
                        currentLink.EncodedBits{iCodeword} = currentLink.ChannelCoder{iCodeword}.encode(currentLink.InputBits{iCodeword});
                    end
                    
                    % If MUST is enabled, get the bits of the farUE.
                    if currentLink.IsMUSTNearUE
                        farUELink = Links{obj.ID, currentLink.MUSTFarUE};
                        farUELink.isScheduled = 1;
                        nCodewordsFar = length(farUELink.layerMappingTable.Downlink{nStreams});
                        farUELink.InputBits = cell(nCodewordsFar,1);
                        farUELink.EncodedBits = cell(nCodewordsFar,1);
                        for iCodeword = 1:nCodewordsFar
                                farUELink.InputBits{iCodeword} = randi([0 1], farUELink.NrInputBits, 1);
                                farUELink.EncodedBits{iCodeword} = farUELink.ChannelCoder{iCodeword}.encode( farUELink.InputBits{iCodeword} );
                        end
                        farUEEncodedBits = farUELink.EncodedBits;
                    else
                        farUEEncodedBits = cell(nCodewords,1);
                    end

                    % Symbol Mapping for this user
                    currentMod  = currentLink.Modulator;
                    spreadSymbolsPerStream = null(1);
                    for iCodeword = 1:nCodewords
                        dataSymbols = currentMod.SignalConstellation.Bit2Symbol(currentLink.EncodedBits{iCodeword}, farUEEncodedBits{iCodeword});
                        currentLink.TransmitSymbols{iCodeword} = dataSymbols;
                        spreadSymbols = currentMod.spreadingTransform(dataSymbols, currentLink.scheduledSubcarriers);
                        spreadSymbolsPerStream = [spreadSymbolsPerStream;reshape(spreadSymbols,currentLink.layerMappingTable.Downlink{nStreams}(iCodeword),[])]; 
                    end
                                        
                    % MIMO precoding of data symbols
                    dataSymbolsPrecoded = currentLink.MIMO.precoding(spreadSymbolsPerStream.', obj.nAntennas);
                    
                    switch currentMod.ChannelEstimationMethod                        
                        case {'PilotAided','Approximate-Perfect','Perfect'}
                            switch currentMod.Waveform
                                case {'OFDM','f-OFDM','WOLA','UFMC','FBMC'}
                                    pilotSymbols = cell(obj.nAntennas, 1);
                                    modulationMat = zeros(currentMod.WaveformObject.Nr.Subcarriers, currentMod.WaveformObject.Nr.MCSymbols, obj.nAntennas);
                                    for iAntenna = 1:obj.nAntennas
                                        pilotSymbolsTemp = currentMod.SignalConstellation.SymbolMapping(randi(currentMod.SignalConstellation.ModulationOrder,[nnz( logical(currentMod.ChannelEstimator.PilotMatrix(:,:,iAntenna))), 1] ));
                                        pilotSymbols{iAntenna} = pilotSymbolsTemp./abs(pilotSymbolsTemp);
                                        
                                        modulationMatTemp   = zeros(currentMod.WaveformObject.Nr.Subcarriers, currentMod.WaveformObject.Nr.MCSymbols);
                                        modulationMatTemp( logical(currentMod.ChannelEstimator.PilotMatrix(:,:,iAntenna)) ) = pilotSymbols{iAntenna};
                                        modulationMatTemp(~logical(sum(currentMod.ChannelEstimator.PilotMatrix,3)))         = dataSymbolsPrecoded(:,iAntenna);
                                        modulationMat(:,:,iAntenna) = modulationMatTemp;
                                    end
                                    currentMod.PilotSymbols = pilotSymbols;
                                case 'GFDM'
                                    % To be added
                                    error('GFDM with channel estimation is not yet implemented!');
                                otherwise
                                    error('Waveform unknown!');
                            end
                        otherwise
                            error('Channel estimation method unkown!');
                    end %switch channel est. method
                    
                    % Both NearUE and FarUE has the same transmit signal
                    if currentLink.IsMUSTNearUE
                        for iCodeword = 1:nCodewords
                            farUELink.TransmitSymbols{iCodeword} = dataSymbols;
                        end
                        farUELink.Modulator.PilotSymbols = pilotSymbols;                       
                    end
                    
                    % modulation for each link
                    transmitSignalUnitPowerPerAntenna = zeros(currentMod.WaveformObject.Nr.SamplesTotal, obj.nAntennas);
                    for iAntenna = 1:obj.nAntennas
                        transmitSignalUnitPowerPerAntenna(:,iAntenna) = currentMod.WaveformObject.Modulation(modulationMat(:,:,iAntenna));
                    end
                    transmitSignalUnitPowerPerUE(:,:,iScheduledUE) = transmitSignalUnitPowerPerAntenna;
                    
                end %is scheduled
                
            end %for iUE
            
            % combine user signals
            transmitSignalUnitPower = sum(transmitSignalUnitPowerPerUE,3);
            
            % rescale average transmit power
            transmitSignal = transmitSignalUnitPower * 10^((obj.txPower-30)/20) * sqrt(totalSubcarriers/totalScheduledSubcarriers);% / sqrt(obj.nAntennas);
            
            % Just copy it to all connected links? Since they all share
            % the same transmit signal (cost: additional memory consumption)
            for iUE = 1:length(obj.ReceiveUE)
                Links{obj.ID, obj.ReceiveUE(iUE)}.applyNonlinearity(transmitSignal);
                Links{obj.ID, obj.ReceiveUE(iUE)}.transmitPower     = obj.txPower;
                Links{obj.ID, obj.ReceiveUE(iUE)}.powerScale        = sqrt(totalSubcarriers/totalScheduledSubcarriers);% / sqrt(obj.nAntennas);
            end             
        end
    end
%% Receive Functions
   methods
        function processReceiveSignal(obj, totalSignal, Links, params)
           % SISO process the users one by one
           for iUE = 1:length(obj.TransmitUE)
                currentLink = Links{obj.TransmitUE(iUE), obj.ID};
                if currentLink.isScheduled && strcmp(currentLink.Type, 'Primary')
                    currentLink.demodulate(totalSignal, params.phy.noisePower);
                    currentLink.decode;
                end
           end
        end 
   end
       
end

