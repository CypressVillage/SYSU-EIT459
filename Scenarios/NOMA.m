% parameter file for the NOMA scenario from "Versatile Mobile Communications Simulation: The Vienna 5G Link Level Simulator"
%
% The setup is as follows: We have two cells (BSs) with 72 subcarriers each.
% The first cell supports two users only and therefore each one gets 36
% subcarriers, while the second cell supports MUST operation and therefore
% can have two groups of two superimposed users (in the power-domain). In 
% each group, one of the superimposed users has a much higher pathloss than 
% the other one (i.e., a cell edge user). The results are obtained over the
% transmit power of the BSs. This scenario then shows the gain/impact of the 
% MUST operation compared to the case when it is not enabled.

% What to expect: UE1, UE2 and will be performing unaffected (orthogonal). 
% UE3 and UE4 will suffer a performance loss at the low SNR regime due to 
% the extra interference from UE5 and UE6, respectively. Overall, cell 2 will 
% perform better since it can squeeze in UE5 and UE6 and therefore will have 
% a higher downlink throughput when the transmit power is sufficiently high.

%% Topology
% Specifiy all the nodes in ascending order with starting
% index of 1 (BS0 or UE0 is not allowed).
scStr.topology.nodes             = ['BS1,BS2,UE1,UE2,UE3,UE4,UE5,UE6'];         

% Primary (desired) links
scStr.topology.primaryLinks      = [ 'BS1:UE1,'...           
                                    'BS1:UE2,'... 
                                    'BS2:UE3,'... 
                                    'BS2:UE4,'... 
                                    'BS2:UE5,'... 
                                    'BS2:UE6'
                                    ];
                                            
% Links for Joint Tranmission and Detection (future work)                         
scStr.topology.jointTxRxLinks               = [''];  

% Interference Links   
scStr.topology.interferenceGeneration       = 'Automatic';
             
scStr.topology.attenuation                  = 300;                      % in dB, set a very high attenuation level to virtually decouple the two BS
scStr.topology.interferingLinks             = [];

%% General Simulation Parameters
% set link types to simulate
scStr.simulation.simulateDownlink           = true;                     % downlink
scStr.simulation.simulateUplink             = false;                    % uplink
scStr.simulation.simulateD2D                = false;                    % device to device links

% Plot options
scStr.simulation.plotResultsFor             = [1];
scStr.simulation.plotOverSNR                = false;
scStr.simulation.saveData                   = false; 

% Define a sweep parameter
scStr.simulation.sweepParam                 = {'simulation.txPowerBaseStation'};  % Define the parameter to sweep over. This can be almost any simulation parameter.
                                                                        % Most likely it will be the pathloss to obtain results over SNR.
                                                                        
scStr.simulation.sweepValue                 = linspace(-30,40,8);       % Define parameter values to sweep over, in dB. A good starting point for the pathloss is 150 to 110
scStr.simulation.applySweepingTo            = [1];

% Number of simulation frames
scStr.simulation.nFrames                    = 100;                      % Number of frames to simulate per sweep value, adjust to obtain sufficiently small confidence intervals.

%% Physical Transmission Parameters
scStr.simulation.centerFrequency            = 2.5e9;                    % center frequency

scStr.simulation.txPowerBaseStation         = [30];                     % per BS; base station total transmit power in dBm
scStr.simulation.txPowerUser                = [30];                     % per UE; user total transmit power in dBm

scStr.simulation.nAntennasBaseStation       = [2];                      % per BS; number of antennas at the base station
scStr.simulation.nAntennasUser              = [2];                      % per UE; number of antennas at the user
scStr.simulation.userVelocity               = [0];                      % per UE; velocity in m/s

% Links to UE1 and UE3 have pathloss of 80, UE2 and UE4 of 90, and UE5 and UE6 of 110 and 115, respectively (cell edge).
scStr.simulation.pathloss                   = [80,90,80,90,110,115];    % per Link, channel pathloss in dB, this is most likely swept over

% Nonlinearity model
scStr.simulation.downlinkNonlinearity       = false;                    % Apply Rapp nonlinear model to downlink channels: either true of false                                                          
scStr.simulation.amplifierOBO               = [1];                      % Amplifier output back-off, per BS, in dB
scStr.simulation.smoothnessFactor           = [3];                      % Smoothness factor for the Rapp model, per BS, >=0

%% Channel Parameters
scStr.channel.dopplerModel                  = 'Discrete-Jakes';        
scStr.channel.timeCorrelation               = false;
scStr.channel.spatialCorrelation            = 'none';            
scStr.channel.nPaths                        = 50;
scStr.channel.powerDelayProfile             = 'PedestrianA';
scStr.channel.K                             = 0;
scStr.channel.delta                         = 1;

%% Channel Estimation and Equalization Parameters
scStr.simulation.channelEstimationMethod    = 'Approximate-Perfect';
scStr.simulation.noisePowerEstimation       = false;                                                          
scStr.simulation.pilotPattern               = 'LTE Downlink';                                                        
scStr.simulation.equalizerType              = 'One-Tap';
scStr.simulation.receiverTypeMIMO           = 'MMSE';

%% MIMO Parameters
% Layer mapping
scStr.layerMapping.mode                     = '5G'; 
scStr.layerMapping.table.Uplink             = {1;2;[1,2]};
scStr.layerMapping.table.Downlink           = {1;2;[1,2]};  

% MIMO mode
scStr.modulation.transmissionMode           = 'CLSM';
scStr.modulation.delayDiversity             = 1;

%% Feedback Parameters
scStr.feedback.delay                        = 0;
scStr.feedback.averager.Type                = 'miesm';
scStr.modulation.cqiTable                   = 0; 

% for the custom transmission mode the following parameters are used to configure the feedback
scStr.feedback.enable                       = false;                    % this parameter is ignored, the feedback is automatically enabled for the CLSM transmission mode
scStr.feedback.pmi                          = false;
scStr.feedback.ri                           = false;
scStr.feedback.cqi                          = true; 
                                                                
scStr.modulation.nStreams                   = [2];                   % per Link; number of active spatial streams
scStr.modulation.mcs                        = [15];                     % parameter is unused
% Per link, precoder selection (used when feedback is disabled)
scStr.modulation.precodingMatrix{1}         = 1/sqrt(2)*ones(2,2);   % per Link; employed precoding matrix

%% Modulation Parameters
% waveform
scStr.modulation.waveform                   = {'OFDM'};
                                                                        
% parameters for FBMC
scStr.modulation.prototypeFilter            = 'PHYDYAS-OQAM';           % unused for OFDM
% Parameters for UFMC  
scStr.modulation.nSubcarriersPerSubband     = [12];                     % number of subcarriers per subband

% time and bandwidth setup (number of subcarriers, frame duration, CP
% length, sampling rate)
scStr.modulation.numerOfSubcarriers         = [72];                     % per BS; number of used subcarriers
scStr.modulation.subcarrierSpacing          = [15e3];                   % per BS; per base station in Hz
scStr.modulation.nSymbolsTotal              = [15];                     % per BS; total number of time-symbols per frame, the frame duration will be nSymbolsTotal/subcarrierSpacing
scStr.modulation.nGuardSymbols              = [1];                      % per BS; select how many of the total time-symbols will be used as guard symbols (cyclic prefix in OFDM)
scStr.modulation.samplingRate               = 'Automatic';              % sampling rate has to be the same for all nodes (across all base stations):
                                                                        % either numeric value for manual setting or 'Automatic'

%% Channel Coding Parameters
% All links are operating with the same coding parameters, enter it only once.
scStr.coding.code                           = {'Turbo'};                                                              
scStr.coding.decoding                       = {'MAX-Log-MAP'};                                                                     
scStr.coding.decodingIterations             = [8];

%% Schedule
% static schedule per base station

% BS1 does Orthogonal Multiple Access
scStr.schedule.fixedScheduleDL{1}           = ['UE1:36,UE2:36'];        % schedule for BS1 Downlink
scStr.schedule.fixedScheduleUL{1}           = [];                       % No uplink for BS1.

% BS2 does MUST operation
scStr.schedule.fixedScheduleDL{2}           = ['UE3:36,UE4:36,UE5:UE3,UE6:UE4'];    % schedule for BS2 Downlink
scStr.schedule.fixedScheduleUL{2}           = [];                                   % No uplink for BS2.
                                                                  
                                                                        
