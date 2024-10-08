% this simulation script simulates FER and throughput curves for CQIs 1 to 15 fot a comparison to the LTE-A downlink link level simulator

close all;
clear;
clc;

cd('..');

% define number of MCS values
nCQI            = 15;
nSNRValues      = 12;
SNRRegionsStart = -10:1.9:18;
SNRRegionDelta  = 4;

%% Setup
% select scenario
simulationScenario = 'LTE_comp_scenario'; 
                                                
% load parameters according to scenario                                               
simParams = Parameters.SimulationParameters( simulationScenario );

% init results
throughputOverall   = zeros(nCQI, 2, nSNRValues);
FEROverall          = zeros(nCQI, 2, nSNRValues);

for cqi = 1:nCQI
    fprintf('------- CQI %d -------\n',cqi);
    
    % set MCS value
    simParams.modulation.mcs = cqi;
    simParams.simulation.sweepValue = 143 - linspace(SNRRegionsStart(cqi), SNRRegionsStart(cqi)+SNRRegionDelta, nSNRValues);

    % generate network topology and links between nodes
    [Links, BS, UE] = Topology.getTopology(simParams);

    nBS         = length(BS);
    nUE         = length(UE);
    dimLinks    = length(Links);
    nFrames     = simParams.simulation.nFrames;

    % initialize result variables
    perSweepResults         = cell(length(simParams.simulation.sweepValue), nFrames);
    simResults              = cell(1, length(simParams.simulation.sweepValue));
    averageFrameDuration    = 0;

    %% Simulation loop
    startTime   = tic;
    myCluster   = parcluster('local');
    NumWorkers  = myCluster.NumWorkers;
    fprintf(['------- Started -------', '\n']);

    % loop over sweep parameter
    parfor iSweep = 1:length(simParams.simulation.sweepValue) % this may be 'for' or 'parfor'
        % update sweep value
        simParams.UpdateSweepValue(iSweep); %#ok

        % based on the updated sweep parameter, regenerate the network
        [Links, BS, UE] = Topology.getTopology(simParams);

        % Objects initialization should be peformed here
        Links = simParams.initializeLinks(Links, BS, UE); %#ok

        % save average frame duration
        if iSweep == 1
            averageFrameDuration(iSweep) = simParams.simulation.averageFrameDuration;
        end

        % Prepare per sweep value results
        perSweepResults = cell(dimLinks, dimLinks, nFrames);

        for iFrame = 1:nFrames
            % link adaptation
            % update all links
            nBS = length( BS );
            L   = length(Links);
            for iBS = 1:nBS
                for iLink = 1:L
                    % Downlink
                    if simParams.simulation.simulateDownlink && ~isempty(Links{BS{iBS}.ID, iLink}) && strcmp(Links{BS{iBS}.ID, iLink}.Type, 'Primary')
                        Links{BS{iBS}.ID, iLink}.updateLink( simParams, Links, iFrame );
                    elseif simParams.simulation.simulateDownlink && ~isempty(Links{BS{iBS}.ID, iLink}) && strcmp(Links{BS{iBS}.ID, iLink}.Type, 'Interference')
                        Links{BS{iBS}.ID, iLink}.Channel.NewRealization(iFrame);
                    end
                    % Uplink
                    if simParams.simulation.simulateUplink && ~isempty(Links{iLink, BS{iBS}.ID}) && strcmp(Links{iLink, BS{iBS}.ID}.Type, 'Primary')
                        Links{iLink, BS{iBS}.ID}.updateLink( simParams, Links,iFrame );
                    elseif simParams.simulation.simulateUplink && ~isempty(Links{iLink, BS{iBS}.ID}) && strcmp(Links{iLink, BS{iBS}.ID}.Type, 'Interference')
                        Links{iLink, BS{iBS}.ID}.Channel.NewRealization(iFrame);
                    end
                end
            end

            %% Downlink
            if simParams.simulation.simulateDownlink
            % All BSs generate their transmit signal for this frame
                for iBS = 1:nBS
                    BS{iBS}.generateTransmitSignal(Links);
                end

                for iUE = 1:nUE
                    UEID = UE{iUE}.ID;
                    primaryLink = Links{UE{iUE}.TransmitBS(1), UEID};             
                    if primaryLink.isScheduled
                        primaryLink.generateReceiveSignal();
                        UETotalSignal = primaryLink.ReceiveSignal;

                        % Collect signals from all other BSs
                        for iBS = 2:length(UE{iUE}.TransmitBS)
                            currentLink = Links{UE{iUE}.TransmitBS(iBS), UEID};
                            currentLink.generateReceiveSignal();
                            UETotalSignal = Channel.addSignals(UETotalSignal, currentLink.ReceiveSignal);
                        end
                        % correct signal length
                        UETotalSignal = Channel.correctSignalLength(UETotalSignal, primaryLink.Modulator.WaveformObject.Nr.SamplesTotal);

                        % add noise
                        UETotalSignal = UETotalSignal + Channel.AWGN( simParams.phy.noisePower, length(UETotalSignal), UE{iUE}.nAntennas );

                        % process received signal
                        UE{iUE}.processReceiveSignal(UETotalSignal, Links, simParams);

                        % Collect the results which is now stored in the primary link
                        primaryLink.calculateSNR(simParams.constants.BOLTZMANN, simParams.phy.temperature);
                        perSweepResults{UE{iUE}.TransmitBS(1), UEID, iFrame} = primaryLink.getResults(simParams.simulation.saveData);
                    end
                end

            end 

            %% Time calculation
            % Some basic time calculation.
            intermediateTime = tic;
            if mod(iFrame, 20) == 0
                fprintf('Sweep: %i/%i, Frame: %i/%i, approx. %.0fs left\n', iSweep,length(simParams.simulation.sweepValue),iFrame,nFrames, double(intermediateTime-startTime)*1e-6*(nFrames*length(simParams.simulation.sweepValue)-(iFrame+(iSweep-1)*nFrames))/(iFrame+(iSweep-1)*nFrames));
            end
        end % for iFrame
        simResults{iSweep} = perSweepResults;
    end % parfor iSweep

    %% post process simulation results
    if simParams.simulation.simulateDownlink
        downlinkResults = Results.SimulationResults( nFrames, length(simParams.simulation.sweepValue), nBS, nUE, averageFrameDuration, 'downlink' );
        downlinkResults.collectResults( simResults, UE );
        downlinkResults.postProcessResults();
        % plot results
        %Results.plotResults( downlinkResults, 'downlink', simParams, UE, BS );
    end

    fprintf(['------- Done -------', '\n']);
    toc(startTime);
    
    % save results for this CQI value
    throughputOverall(cqi,1,:)  = mean(downlinkResults.userResults(1).SNR);
    throughputOverall(cqi,2,:)  = downlinkResults.userResults(1).throughput.mean;
    
    FEROverall(cqi,1,:)         = mean(downlinkResults.userResults(1).SNR);
    FEROverall(cqi,2,:)         = downlinkResults.userResults(1).FER.mean;
    
end

save('results/5GSimReferenceCurves','throughputOverall','FEROverall');
