% parameter file for the multi-link simulation scenario from "Versatile Mobile Communications Simulation: The Vienna 5G Link Level Simulator"
%
% you may consider to change the waveform parameter to generate all result curves for this scenario

%% Topology
scStr.topology.nodes                        = ['BS1,BS2,UE1,UE2'];      % Specify all nodes in the network

% Primary (desired) links
scStr.topology.primaryLinks                 = ['UE1:BS1,' ...           % two cells with single user each
                                               'UE2:BS2'];              
                                                
% Links for Joint Tranmission and Detection (future work)                         
scStr.topology.jointTxRxLinks               = [''];  

% Interference Links   
scStr.topology.interferenceGeneration       = 'Automatic';              % 'Automatic' automatically generates all possible interference links
                                                                        % in this scenario the interfering links are UE1:BS2 and UE2:BS1
                                                                        
scStr.topology.attenuation                  = 107;                      % the interfering link attenuation is set to the same value as the desired link's path loss
scStr.topology.interferingLinks             = [];                       % not neccessary for automatic interference link generation
                                                                        
%% General Simulation Parameters
% Set link types to simulate
scStr.simulation.simulateDownlink           = false;
scStr.simulation.simulateUplink             = true;                     % Uplink only
scStr.simulation.simulateD2D                = false;

% Plot options
scStr.simulation.plotResultsFor             = [1];
                                                                         
scStr.simulation.plotOverSNR                = true;                                                                 
scStr.simulation.saveData                   = false;

% Define a sweep parameter
scStr.simulation.sweepParam                 = {'simulation.txPowerUser'};  % sweep over a user's transmit power
                                                                                                                                        
scStr.simulation.sweepValue                 = 10:5:60;                  % transmit power of UE2

scStr.simulation.applySweepingTo            = [0,1];                    % apply sweep to second user only
% Number of simulation frames
scStr.simulation.nFrames                    = 2000;                     % Number of frames to simulate per sweep value, adjust to obtain sufficiently small confidence intervals.

%% Physical Transmission Parameters
scStr.simulation.centerFrequency            = 2.5e9;                    % center frequency

scStr.simulation.txPowerBaseStation         = [30,30];                  % per BS; BS total transmit power in dBm
scStr.simulation.txPowerUser                = [30,30];                  % per UE; UE total transmit power in dBm

scStr.simulation.nAntennasBaseStation       = [1];                      % singel antenna base stations
scStr.simulation.nAntennasUser              = [1];                      % single antenna users
scStr.simulation.userVelocity               = [0];                      % no user velocity, static scenario

scStr.simulation.pathloss                   = [107];                    % path loss in dB for an SNR of approx. 40dB

% Nonlinearity model
scStr.simulation.downlinkNonlinearity       = false;                   % Apply Rapp nonlinear model to downlink channels    {'Yes'/'No'}                                                         
scStr.simulation.amplifierOBO               = [1];                      % Amplifier output back-off, per BS, in dB
scStr.simulation.smoothnessFactor           = [3];                      % Smoothness factor for the Rapp model, per BS, >=0

%% Channel Parameters
scStr.channel.dopplerModel                  = 'Discrete-Jakes';         % the Doppler model is not of interest as the user velocity is zero and the scenario is static                                                                        
scStr.channel.timeCorrelation               = false;                    % no time correlaltion                        
scStr.channel.spatialCorrelation            = 'none';                   % no spatial correlation                                                              
scStr.channel.nPaths                        = 50;                       % number of propagation paths for the Doppler model
                                                                        % for time correlated fading 10 is sufficient, otherwise 50 or higher

scStr.channel.powerDelayProfile             = 'PedestrianA';            % a channel model with a low RMS delay spread is chosen such that there is no ISI or ICI

scStr.channel.K                             = 0;                       	% Rayleigh fading
scStr.channel.delta                         = 1;

%% Channel Estimation and Equalization Parameters
scStr.simulation.channelEstimationMethod    = 'Approximate-Perfect';    % assume perfect channel knowledge for this scenario
scStr.simulation.noisePowerEstimation       = false;                         
scStr.simulation.pilotPattern               = 'Diamond';             
scStr.simulation.equalizerType              = 'One-Tap';             
scStr.simulation.receiverTypeMIMO           = 'ZF';

%% MIMO Parameters
% Layer mapping
scStr.layerMapping.mode                     = '5G';                     % LTE, 5G or custom
scStr.layerMapping.table.Uplink             = {1;2;[1,2]};
scStr.layerMapping.table.Downlink           = {1;2;[1,2]};

% MIMO mode
scStr.modulation.transmissionMode           = 'custom';                 % this mode allows manual settings for the feedback calculation                                                                  
scStr.modulation.delayDiversity             = 1;                        % 0 no delay 1 large delay diversity (only for OLSM) 

%% Feedback Parameters
% feedback is disabled for this scenario
scStr.feedback.delay                        = 0;
scStr.feedback.averager.Type                = 'miesm';
scStr.modulation.cqiTable                   = 0;                        % no feedback is used for this sceario
scStr.feedback.enable                       = false; 
scStr.feedback.pmi                          = false;                        
scStr.feedback.ri                           = false;                       
scStr.feedback.cqi                          = false;  

% RI, PMI and CQI setting
scStr.modulation.nStreams                   = [1, 1];                   % single antenna, SISO transmission
scStr.modulation.precodingMatrix{1}         = 1;                        % no spatial precoding matrix required
scStr.modulation.precodingMatrix{2}         = 1;
scStr.modulation.mcs                        = [12,12];                  % per link.  

%% Modulation Parameters
% Waveform
scStr.modulation.waveform                   = {'OFDM'};                 % you may want to try 'f-OFDM' or 'FBMC' in this scenario
                                                                        
% Parameters for FBMC
scStr.modulation.prototypeFilter            = 'PHYDYAS-OQAM';           % prototype filter for FBMC
% Parameters for UFMC  
scStr.modulation.nSubcarriersPerSubband     = [12];                     % number of subcarriers per subband
                                                                 
% Time and bandwidth setup (number of subcarriers, frame duration, CP length, sampling rate)
scStr.modulation.numerOfSubcarriers         = [72,36];                  % per BS; number of used subcarriers
scStr.modulation.subcarrierSpacing          = [15e3,30e3];              % per BS; per base station in Hz
scStr.modulation.nSymbolsTotal              = [15,30];                  % per BS; total number of time-symbols per frame, the frame duration will be nSymbolsTotal/subcarrierSpacing
scStr.modulation.nGuardSymbols              = [1,2];                    % per BS; select how many of the total time-symbols will be used as guard symbols (cyclic prefix in OFDM)
scStr.modulation.samplingRate               = 'Automatic';                % sampling rate has to be the same for all nodes (across all base stations):
                                                                        % either numeric value for manual setting or 'Automatic'

%% Channel Coding Parameters
scStr.coding.code                           = {'LDPC'};                 % LDPC channel coding
scStr.coding.decoding                       = {'PWL-Min-Sum'};                            
scStr.coding.decodingIterations             = [8];                      % number of decoding iterations

%% Schedule
% static schedule per base station
scStr.schedule.fixedScheduleDL{1}           = [];                       % schedule for BS1 downlink
scStr.schedule.fixedScheduleUL{1}           = ['UE1:34,none:38'];       % schedule for BS1 uplink
scStr.schedule.fixedScheduleUL{2}           = ['none:19,UE2:17'];       % schedule for BS2 uplink
                    

