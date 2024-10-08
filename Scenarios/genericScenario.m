% parameter file for the generic scenario
%
% The idea behind this scenario is, to show all possible options for all
% parameters. It is also a good starting point to define your own
% simulation scenarios.

% You can either specify the parameter for each node/link individually or if 
% you want all of them to have the same parameter, then you just enter a
% single value. E.g. scStr.simulation.userVelocity = [0]; means that all
% users will have the same velocity of zero m/s.

%% Topology
% Specifiy all the nodes in ascending order with starting index of 1 (BS0 or UE0 is not allowed).
scStr.topology.nodes                        = ['BS1,BS2,UE1,UE2'];          % Specify all nodes in the network

% Primary (desired) links
scStr.topology.primaryLinks                 = ['BS1:UE1,' ...
                                               'BS2:UE2' ...
                                               ];              % Links to be considered as disired links
                                                
                                            
% Links for Joint Tranmission and Detection (future work)                         
scStr.topology.jointTxRxLinks               = [''];  

% Interference Links   
scStr.topology.interferenceGeneration       = 'Automatic';              % Generation of interference links
                                                                        % 'Automatic' automatically generates all possible interference links
                                                                        % 'Manual' manually choose interference links and their attenuation
                                                            
scStr.topology.attenuation                  = 200;                      % In dB, for all links under 'Automatic' generation

scStr.topology.interferingLinks             = [ ];                      % For 'Manual' generation, specifiy the interference links
                                                                        % Value after the * symbol denotes the attenuation in dB.  
                                                                        % BS1:UE3*20 means that the link from BS1 to UE3 is an
                                                                        % interference link with attenuation equal to 20 dB;
                                                                        
%% General Simulation Parameters
% Set link types to simulate
scStr.simulation.simulateDownlink           = true;                     % Downlink
scStr.simulation.simulateUplink             = false;                    % Uplink
scStr.simulation.simulateD2D                = false;                    % Device-To-Device (D2D) links

% Plot options
scStr.simulation.plotResultsFor             = [1];                      % Choose the nodes that will get their results plotted (1 or 0 for each node).
                                                                        % Enter a single 1 if you want the results of all nodes to be plotted, or 0 if you want no plots.
                                                                    
                                                                        
                                                                         
scStr.simulation.plotOverSNR                = true;                     % Select to plot results over SNR instead of pathloss (if pathloss is selected as the sweep parameter).                                                                        
scStr.simulation.saveData                   = false;                    % Set to true to save bits and symbols.

% Define a sweep parameter
scStr.simulation.sweepParam                 = {'simulation.pathloss'};  % Define the parameter to sweep over. This can be almost any simulation parameter.
                                                                        % Most likely it will be the pathloss to obtain results over SNR.
                                                                                                                                        
scStr.simulation.sweepValue                 = linspace(150,90,6);       % Define parameter values to sweep over. A good starting point for the pathloss is 150 to 110.

scStr.simulation.applySweepingTo            = [1];                      % Define the nodes (or links, depending on the sweep parameter) on which the sweep parameter will be applied (1 or 0 for each node/link).
                                                                        % All other nodes (or links) set to 0 here will use the custom values entered below. 
                                                                        % Enter a single 1 if you want all of them to be swept over.
% Number of simulation frames
scStr.simulation.nFrames                    = 100;                      % Number of frames to simulate per sweep value, adjust to obtain sufficiently small confidence intervals.

%% Physical Transmission Parameters
scStr.simulation.centerFrequency            = 2.5e9;                    % center frequency

scStr.simulation.txPowerBaseStation         = [30];                     % per BS; BS total transmit power in dBm
scStr.simulation.txPowerUser                = [30,30];                  % per UE; UE total transmit power in dBm

scStr.simulation.nAntennasBaseStation       = [2];                      % per BS; number of antennas at the BS
scStr.simulation.nAntennasUser              = [1];                      % per UE; number of antennas at the BS
scStr.simulation.userVelocity               = [0];                      % per UE; velocity in m/s

scStr.simulation.pathloss                   = [80];                     % per link, channel pathloss in dB, this is most likely swept over

% Nonlinearity model
scStr.simulation.downlinkNonlinearity       = false;                    % Apply Rapp nonlinear model to downlink channels: either true of false       
scStr.simulation.amplifierOBO               = [1];                      % Amplifier output back-off, per BS, in dB
scStr.simulation.smoothnessFactor           = [3];                      % Smoothness factor for the Rapp model, per BS, >=0

%% Channel Parameters
scStr.channel.dopplerModel                  = 'Discrete-Jakes';         % implementation of a Doppler spectrum:
                                                                        % 'Jakes'
                                                                        % 'Uniform'
                                                                        % 'Discrete-Jakes'
                                                                        % 'Discrete-Uniform'
                                                                        
scStr.channel.timeCorrelation               = false;                    % channel time correlation between frames
                           
scStr.channel.spatialCorrelation            = 'none';                   % spatial correlation according to TS36.101, Annex B
                                                                        % 'none': spatial correlation model is deactivated
                                                                        % 'low': no spatial correlation
                                                                        % 'medium': correlation at UE and BS, but higher at UE
                                                                        % 'high': high correlation at UE and at BS
                                                                
scStr.channel.nPaths                        = 50;                       % number of propagation paths for the Doppler model
                                                                        % for time correlated fading 10 is sufficient, otherwise 50 or higher

scStr.channel.powerDelayProfile             = 'PedestrianA';            % power delay profile model, possible choices are:
                                                                        % 'AWGN'
                                                                        % 'Flat'
                                                                        % 'PedestrianA'
                                                                        % 'PedestrianB'
                                                                        % 'VehicularA'
                                                                        % 'VehicularB'
                                                                        % 'ExtendedPedestrianA'
                                                                        % 'ExtendedVehicularA'
 
                                                                       
% Fading parameters for the time-invariant case for all users: TWDP model
% These parameters will be carried out only if all users are stationary, i.e the userVelocity parameter is zero for all users
% If there is a mix of moving and stationary users, Rayleigh fading will be  assigned to the stationary ones!
scStr.channel.K                             = 0;                        % TWDP model parameter K,      K > 0
scStr.channel.delta                         = 1;                        % TWDP model parameter delta,  0 >= delta >= 1
                                                                        % TWDP parameters and the represented statistic:
                                                                        % K >>  ; delta = 0       : No fading
                                                                        % K > 0 ; delta = 0       : Rician fading
                                                                        % K = 0 ; delta arbitrary : Rayleigh fading
                                                                        % K >>  ; delta ~ 1       : Hyper-rayleigh fading 

%% Channel Estimation and Equalization Parameters
scStr.simulation.channelEstimationMethod    = 'Approximate-Perfect';    % channel estimation method: 
                                                                        % 'Approximate-Perfect'
                                                                        % 'PilotAided'
                                                                        
scStr.simulation.noisePowerEstimation       = false;                    % if the noise and interference power are estimated
                                                                        % true: estimate interference and noise power at RX
                                                                        % false: use perfect knowledge of noise power, neglect interference power
                                                                        
scStr.simulation.pilotPattern               = 'LTE Downlink';           % pilot symbol allocation pattern
                                                                        % 'Rectangular'
                                                                        % 'Diamond'
                                                                        % 'LTE Downlink'
                                                                    
scStr.simulation.equalizerType              = 'One-Tap';                % Multicarrier equalizer: 'One-Tap'

scStr.simulation.receiverTypeMIMO           = 'ZF';                     % MIMO detector: 'ZF','MMSE','Sphere','ML'

%% MIMO Parameters
% Layer mapping
scStr.layerMapping.mode                     = 'LTE';                     % 'LTE', '5G' or 'custom'

scStr.layerMapping.table.Uplink             = {1;2;[1,2]};
scStr.layerMapping.table.Downlink           = {1;2;[1,2]};              % table{n} represents the codeword to layer mapping for n layers,
                                                                        % table{3} = [1,2] the number of enteries in the array equals the number of
                                                                        % codewords and each entery represents the number of layers the codeword it
                                                                        % mapped to.
% MIMO mode
scStr.modulation.transmissionMode           = 'custom';                 % 'OLSM' open loop spatial multiplexing
                                                                        % 'CLSM' closed loop spatial multiplexing
                                                                        % 'TxD' transmit diversity
                                                                        % 'custom' for the custom mode the feedback has to be configured manually (see scStr.feedback)
                                                                                                                                           
scStr.modulation.delayDiversity             = 1;                        % 0 no delay 1 large delay diversity (only for OLSM) 

%% Feedback Parameters
scStr.feedback.delay                        = 0;                        % Feedback Delay

scStr.feedback.averager.Type            	= 'miesm';                  % Method used to calculate an effective SNR
                                                                        % 'eesm'
                                                                        % 'miesm'
                                                                        
scStr.modulation.cqiTable                   = 0;                        % CQI table selection
                                                                        % 0: TS 38.214 Table 5.2.2.1-2 (QPSK to 64 QAM)
                                                                        % 1: TS 38.214 Table 5.2.2.1-3 (QPSK to 256 QAM)
                                                                        
% for the custom transmission mode the following parameters are used to configure the feedback
scStr.feedback.enable                       = false;                    % set this to ture to enable the feedback
% when the feedback is enabled the following parameters are used to configure the individual indicators:
scStr.feedback.pmi                          = false;
scStr.feedback.ri                           = false; 
scStr.feedback.cqi                          = true;

% when feedback is disabled, the following parameters are used:                                                                         
scStr.modulation.nStreams                   = [1, 1];                   % per Link; number of active spatial streams

% Per link, precoder selection. Precoding matrix needs to be of size [nStreams x nAntennas]
scStr.modulation.precodingMatrix{1}         = 1/sqrt(2)*ones(2,1);      % Link 1 employed precoding matrix
scStr.modulation.precodingMatrix{2}         = 1/sqrt(2)*ones(2,1);      % Link 2 employed precoding matrix

% CQI selection
scStr.modulation.mcs                        = [8,8];                  % 1 to 15, per link


%% Modulation Parameters
% Waveform
scStr.modulation.waveform                   = {'f-OFDM','OFDM'};        % per BS; employed transmission waveform:
                                                                        % 'OFDM' - Orthogonal Frequency Division Multiplexing
                                                                        % 'FBMC' - Filer Bank MultiCarrier
                                                                        % 'UFMC' - Universal Filtered MultiCarrier
                                                                        % 'f-OFDM' - filtered OFDM
                                                                        % 'WOLA' - Weighted OverLap and Add
                                                                        
% Parameters for FBMC
scStr.modulation.prototypeFilter            = 'PHYDYAS-OQAM';           % prototype filter for FBMC:
                                                                        % 'Hermite-OQAM'
                                                                        % 'Hermite-QAM'
                                                                        % 'Rectangle-QAM'
                                                                        % 'RRC-OQAM'
                                                                        % 'RRC-QAM'
                                                                        % 'PHYDYAS-OQAM'
                                                                        % 'PHYDYAS-QAM'
% Parameters for UFMC  
scStr.modulation.nSubcarriersPerSubband     = [12];                     % number of subcarriers per subband
                                                                 
% Numerology setup
scStr.modulation.numerOfSubcarriers         = [72,72];                  % per BS; number of used subcarriers
scStr.modulation.subcarrierSpacing          = [30e3,15e3];              % per BS; per base station in Hz
scStr.modulation.nSymbolsTotal              = [30,15];                  % per BS; total number of time-symbols per frame, the frame duration will be nSymbolsTotal/subcarrierSpacing
scStr.modulation.nGuardSymbols              = [2,1];                    % per BS; select how many of the total time-symbols will be used as guard symbols (cyclic prefix in OFDM)
scStr.modulation.samplingRate               = 'Automatic';              % ther is a single sampling rate for all nodes (all users and all base stations):
                                                                        % either numeric value for manual setting or 'Automatic'

%% Channel Coding Parameters
scStr.coding.code                           = {'Turbo'};                % per Link; channel code:
                                                                        % 'Turbo'
                                                                        % 'LDPC'
                                                                        % 'TB-Convolutional'
                                                                        % 'Polar'
                                                                    
scStr.coding.decoding                       = {'Linear-Log-MAP'};       % per Link; decoding algorithm:
                                                                        % TB-Convolutional: 'Log-MAP', 'MAX-Log-MAP'
                                                                        % Turbo: 'Log-MAP', 'MAX-Log-MAP', 'Linear-Log-MAP'
                                                                        % LDPC: 'Sum-Product', 'Min-Sum', 'PWL-Min-Sum'
                                                                        % Polar: 'SC', 'List-SC', 'CRC-List-SC'
                                                                            
scStr.coding.decodingIterations             = [8,8];                    % per Link; number of decoding iterations or list size (typically 8 for turbo, 16 for LDPC, 8 for polar)

%% Schedule
% static schedule per base station
scStr.schedule.fixedScheduleDL{1}           = ['UE1:72'];               % schedule for BS1 downlink
scStr.schedule.fixedScheduleDL{2}           = ['UE2:72'];               % schedule for BS1 downlink
scStr.schedule.fixedScheduleUL{1}           = ['UE1:36,none:36'];       % schedule for BS1 uplink
scStr.schedule.fixedScheduleUL{2}           = ['UE2:36,none:36'];       % schedule for BS2 uplink

                                                                        



