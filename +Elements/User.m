classdef User < handle & Elements.Node
%% The User Class (User)
%   Authors
%     - ...
%     - ...
%   (c) 2016 Institute of Telecommunications, TU Wien.
%   www.nt.tuwien.ac.at 
%% User Properties
    properties
        TransmitSignal
        txPower
        nAntennas
    end
%% The User Functions
%% Initialization
    methods
        function obj = User(name, id, nAntennas, transmitPower, plotResults)
            obj = obj@Elements.Node(name, id, 'User', plotResults);
            obj.nAntennas   = nAntennas;
            obj.txPower     = transmitPower;
        end
    end
%% Transmit Functions
    methods
        function generateTransmitSignal(obj, Links)
            % Uplink from the UE to its primary BS
            primaryLink = Links{obj.ID, obj.ReceiveBS(1)};
            if any(primaryLink.schedule(:))
                primaryLink.isScheduled = 1;
                nStreams = primaryLink.MIMO.nStreams;
                nCodewords = length(primaryLink.layerMappingTable.Uplink{nStreams});
                primaryLink.InputBits   = cell(nCodewords,1);
                codedBits = cell(nCodewords,1);
                dataSymbols = cell(nCodewords,1);
                spreadSymbolsPerStream = null(1);
                
                primaryMod = primaryLink.Modulator;

                for iCodeword = 1:nCodewords
                    primaryLink.InputBits{iCodeword} =  randi([0 1], primaryLink.NrInputBits(iCodeword), 1);
                    % Channel Coding
                    codedBits{iCodeword} = primaryLink.ChannelCoder{iCodeword}.encode(primaryLink.InputBits{iCodeword});
                    % Symbol Mapping
                    dataSymbols{iCodeword} = primaryMod.SignalConstellation.Bit2Symbol(codedBits{iCodeword});
                    % Data symbol spreading
                    spreadSymbols = primaryMod.spreadingTransform(dataSymbols{iCodeword}, primaryLink.scheduledSubcarriers);
                    % Layer mapping
                    spreadSymbolsPerStream = [spreadSymbolsPerStream;reshape(spreadSymbols,primaryLink.layerMappingTable.Uplink{nStreams}(iCodeword),[])];                    
                    
                end

                 
                
                % MIMO precoding
                dataSymbolsPrecoded = primaryLink.MIMO.precoding( spreadSymbolsPerStream.', obj.nAntennas );
                
                switch primaryMod.ChannelEstimationMethod
                    case {'PilotAided','Approximate-Perfect','Perfect'}
                        switch primaryMod.Waveform
                            case {'OFDM','f-OFDM','WOLA','UFMC','FBMC'}
                                modulationMat = zeros(primaryLink.Modulator.WaveformObject.Nr.Subcarriers, primaryLink.Modulator.WaveformObject.Nr.MCSymbols, obj.nAntennas);
                                pilotSymbols = cell(obj.nAntennas, 1);
                                for iAntenna = 1:obj.nAntennas
                                    pilotSymbolsTemp = primaryMod.SignalConstellation.SymbolMapping(randi(primaryMod.SignalConstellation.ModulationOrder,[nnz( logical(primaryMod.ChannelEstimator.PilotMatrix(:,:,iAntenna)) ) 1]));
                                    pilotSymbols{iAntenna} = pilotSymbolsTemp./abs(pilotSymbolsTemp);
                                    
                                    modulationMatTemp = zeros(primaryLink.Modulator.WaveformObject.Nr.Subcarriers, primaryLink.Modulator.WaveformObject.Nr.MCSymbols);
                                    modulationMatTemp(logical(primaryMod.ChannelEstimator.PilotMatrix(:,:,iAntenna)))   = pilotSymbols{iAntenna};
                                    modulationMatTemp(~logical(sum(primaryMod.ChannelEstimator.PilotMatrix,3)))         = dataSymbolsPrecoded(:,iAntenna);
                                    modulationMat(:,:,iAntenna) = modulationMatTemp;
                                end
                                primaryMod.PilotSymbols = pilotSymbols;
%                             case 'FBMC'
%                                 % To be added
%                                 error('FBMC with channel estimation is not yet implemented!');
                            case 'GFDM'
                                % To be added
                                error('GFDM with channel estimation is not yet implemented!');
                            otherwise
                                error('Waveform unknown!');
                        end
                    otherwise
                        error('Channel estimation method unkown!');
                end
                
                transmitSignalUnitPower = zeros( primaryMod.WaveformObject.Nr.SamplesTotal, obj.nAntennas ) ;
                for iAntenna = 1:obj.nAntennas
                    transmitSignalUnitPower(:,iAntenna) = primaryMod.WaveformObject.Modulation(modulationMat(:,:,iAntenna));
                end
                
                % normalize transmit signal to tranmit power
                transmitSignal = transmitSignalUnitPower * 10^((obj.txPower-30)/20) * 1/sqrt(obj.nAntennas);
            
                % Just copy it to all connected links? Since they all share
                % the same transmit signal (cost: additional memory consumption)
                for iLink = 1:length(Links)
                    currentLink = Links{obj.ID, iLink};
                    if ~isempty(currentLink)
                        currentLink.isScheduled     = 1;                % For uplink interefence
                        currentLink.EncodedBits     = codedBits;
                        currentLink.TransmitSymbols = dataSymbols;
                        currentLink.TransmitSignal  = transmitSignal;
                        currentLink.powerScale      = 1/sqrt(obj.nAntennas);
                        currentLink.transmitPower   = obj.txPower;
                    end
                end  
            end
        end 
    end
%% Receive Functions
    methods
        function processReceiveSignal(obj, totalSignal, Links, params)
            % Direct processing using the user primary link
            Links{obj.TransmitBS(1), obj.ID}.demodulate(totalSignal, params.phy.noisePower );
            Links{obj.TransmitBS(1), obj.ID}.decode;
        end
    end  
end

