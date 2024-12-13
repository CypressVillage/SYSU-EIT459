%% The Vienna 5G Link Level Simulator v1.1
% www.tc.tuwien.ac.at/vccs
% Please refer to the user manual to get familiar with the simulator structure and mechanics. 
% For any questions, consider our forum www.tc.tuwien.ac.at/forum.

% Main simulator script

close all;
clear;
clc;

%% Setup
% select scenario
simulationScenario = 'LTEAcompliant';           % select a simulation scenario:
                                                % 'genericScenario'
                                                % 'LTEAcompliant'
                                                % 'multiLink'
                                                % 'flexibleNumerology'
                                                % 'NOMA'
                                                                             
% load parameters according to scenario                                               
simParams = Parameters.SimulationParameters( simulationScenario );

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
for iSweep = 1:length(simParams.simulation.sweepValue) % this may be 'for' or 'parfor'
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
                    
                    % % 4.2 用保存的信号文件作为接收机输入
                    % UETotalSignal = primaryLink.TransmitSignal;

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
                    Links{UE{iUE}.TransmitBS(1), UEID}.UETotalSignal_ = UETotalSignal;

                    % Collect the results which is now stored in the primary link
                    primaryLink.calculateSNR(simParams.constants.BOLTZMANN, simParams.phy.temperature);
                    perSweepResults{UE{iUE}.TransmitBS(1), UEID, iFrame} = primaryLink.getResults(simParams.simulation.saveData);
                end
            end
            
         end 
         %% Uplink
         if simParams.simulation.simulateUplink
            % All UEs generate their transmit signal for this frame
            for iUE = 1:nUE
                UE{iUE}.generateTransmitSignal(Links);
            end
            for iBS = 1:nBS
                BSID = BS{iBS}.ID;
                BSTotalSignal = [];
                % Collect signals from all users transmitting to this BS
                % (both primary and interfering users)
                for iUE = 1:length(BS{iBS}.TransmitUE)
                    currentLink = Links{BS{iBS}.TransmitUE(iUE), BSID};
                    if strcmp(currentLink.Type, 'Primary')
                        signalLength = currentLink.Modulator.WaveformObject.Nr.SamplesTotal;
                    end
                    if currentLink.isScheduled
                        currentLink.generateReceiveSignal();                      
                        BSTotalSignal = Channel.addSignals(BSTotalSignal, currentLink.ReceiveSignal);
                    end
                end
                % correct signal length
                BSTotalSignal = Channel.correctSignalLength(BSTotalSignal, signalLength);      
                
                % add noise
                BSTotalSignal = BSTotalSignal + Channel.AWGN( simParams.phy.noisePower, length(BSTotalSignal), BS{iBS}.nAntennas );
                
                % process received signal
                BS{iBS}.processReceiveSignal(BSTotalSignal, Links, simParams);

                % Collect the results
                for iUE = 1:length(BS{iBS}.TransmitUE)
                    currentLink = Links{BS{iBS}.TransmitUE(iUE), BSID};
                    if currentLink.isScheduled && strcmp(currentLink.Type, 'Primary')
                        currentLink.calculateSNR(simParams.constants.BOLTZMANN, simParams.phy.temperature);
                        perSweepResults{BS{iBS}.TransmitUE(iUE), BSID, iFrame} = currentLink.getResults(simParams.simulation.saveData);
                    end
                end  
            end
        end
        %% Device-to-Device
%         if simParams.simulation.simulateD2D
%            % Not yet implemented, but it can be done in a similar fashion
%         end

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
    if sum(simParams.simulation.plotResultsFor) ~= 0
        % Results.plotResults( downlinkResults, 'downlink', simParams, UE, BS );
    end
end
if simParams.simulation.simulateUplink
    uplinkResults = Results.SimulationResults( nFrames, length(simParams.simulation.sweepValue), nBS, nUE, averageFrameDuration, 'uplink' );
    uplinkResults.collectResults( simResults, UE );
    uplinkResults.postProcessResults();
    % plot results
    if sum(simParams.simulation.plotResultsFor) ~= 0
        Results.plotResults( uplinkResults, 'uplink', simParams, UE, BS );
    end
end

%% save results
% generate timestamp
tmpStr = datestr(now);
tmpStr = strrep(tmpStr,':','_');
timeStamp = strrep(tmpStr,' ','_');
% save results (complete workspace)
% save(['./results/results_',timeStamp]);

fprintf(['------- Done -------', '\n']);
toc(startTime);


% % 保存发射信号
% var4_1 = Links{1,2}.TransmitSignal(:, 1);
% save('TransmitSignal.mat', 'var4_1')

%% 保存吞吐量、误码率到.txt文件
% % 保存吞吐量
% throughput = mean(downlinkResults.userResults.throughput.values);
% % 保存误码率
% BER = mean(downlinkResults.userResults.BERUncoded.values);
% fprintf(2, 'Throughput: %f\n', throughput);
% fprintf(2, 'BER: %f\n', BER);
% return

%% 绘图
close all;
for iBS = 1:nBS
    BSID = BS{iBS}.ID;
    for iUE = 1:nUE
        UEID = UE{iUE}.ID;

        if isempty(Links{BSID, UEID}.Modulator)
            continue
        end
        % 绘制接收信号星座图
        figure(iBS*1000+iUE*100+1)
        set(gcf, 'Name', '接收信号星座图');
        scatter(real(Links{BSID, UEID}.Modulator.rxData_(:, 1)), imag(Links{BSID, UEID}.Modulator.rxData_(:, 1)), '.')
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' received signal constellation'])

        % 绘制误码率
        figure(iBS*1000+iUE*100+2)
        set(gcf, 'Name', '误码率')
        bar(downlinkResults.userResults(iUE).BERUncoded.values)
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' Coded BER'])
        xlabel('frame number')
        ylabel('BER Uncoded')

        % 绘制吞吐量
        figure(iBS*1000+iUE*100+3)
        set(gcf, 'Name', '吞吐量')
        bar(downlinkResults.userResults(iUE).throughput.values)
        xlabel('frame number')
        ylabel('throughput')
        title(['BS ', num2str(iBS), ' User ', num2str(iUE)])

        % 绘制发送信号功率谱
        figure(iBS*1000+iUE*100+4)
        set(gcf, 'Name', '发送信号功率谱')
        [pxx, f] = pwelch(Links{BSID, UEID}.TransmitSignal(:, 1), [], [], [], simParams.modulation.samplingRate);
        plot(f, 10*log10(pxx))
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' transmit signal power spectrum'])
        xlabel('frequency/Hz')
        ylabel('power/dB')

        % 绘制接收信号功率谱
        figure(iBS*1000+iUE*100+5)
        set(gcf, 'Name', '接收信号功率谱')
        [pxx, f] = pwelch(Links{BSID, UEID}.UETotalSignal_(:, 1), [], [], [], simParams.modulation.samplingRate);
        plot(f, 10*log10(pxx))
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' received signal power spectrum'])
        xlabel('frequency/Hz')
        ylabel('power/dB')

        % 绘制发送信号时域波形
        figure(iBS*1000+iUE*100+6)
        set(gcf, 'Name', '发送信号时域波形')
        t = 0:1/simParams.modulation.samplingRate:(length(Links{BSID, UEID}.TransmitSignal(:, 1))-1)/simParams.modulation.samplingRate;
        subplot(2,1,1)
        plot(t, abs(Links{BSID, UEID}.TransmitSignal(:, 1)))
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' transmit signal amplitude'])
        xlabel('time/s')
        ylabel('amplitude')
        subplot(2,1,2)
        plot(t, angle(Links{BSID, UEID}.TransmitSignal(:, 1)))
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' transmit signal phase'])
        xlabel('time/s')
        ylabel('phase')

        % 绘制接收信号时域波形
        figure(iBS*1000+iUE*100+7)
        set(gcf, 'Name', '接收信号时域波形')
        t = 0:1/simParams.modulation.samplingRate:(length(Links{BSID, UEID}.UETotalSignal_(:, 1))-1)/simParams.modulation.samplingRate;
        subplot(2,1,1)
        plot(t, abs(Links{BSID, UEID}.UETotalSignal_(:, 1)))
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' received signal amplitude'])
        xlabel('time/s')
        ylabel('amplitude')
        subplot(2,1,2)
        plot(t, angle(Links{BSID, UEID}.UETotalSignal_(:, 1)))
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' received signal phase'])
        xlabel('time/s')
        ylabel('phase')

        % 绘制真实信道
        figure(iBS*1000+iUE*100+8)
        set(gcf, 'Name', '真实信道')
        channel = Links{BSID, UEID}.Modulator.Channel(:,:,1);
        x = 1:1:size(channel, 1);
        y = 1:1:size(channel, 2);
        [X, Y] = meshgrid(x, y);
        surf(X, Y, 10*log(abs(channel')))
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' real channel'])
        xlabel('subcarrier')
        ylabel('OFDM symbol')
        zlabel('channel gain')

        % 绘制估计信道
        figure(iBS*1000+iUE*100+9)
        set(gcf, 'Name', '估计信道')
        channel = Links{BSID, UEID}.Modulator.perfectChannel_(:,:,1);
        surf(X, Y, 10*log(abs(channel')))
        title(['BS ', num2str(iBS), ' User ', num2str(iUE), ' estimated channel'])
        xlabel('subcarrier')
        ylabel('OFDM symbol')
        zlabel('channel gain')

        Links{BSID, UEID}.Channel.PlotTimeCorrelation(Links{BSID, UEID}.Modulator.WaveformObject.Implementation.TimeSpacing);
        [FCF_real, frequency] = Links{BSID, UEID}.Channel.PlotFrequencyCorrelation(2);
        [PDP_real, tau] = Links{BSID, UEID}.Channel.PlotPowerDelayProfile();
        
        CFR = Links{BSID, UEID}.Modulator.Channel(:,:,1);
        CIR = ifft(CFR);

        PDP_estimated = mean(abs(CIR).^2, 2);
        PDP_estimated = PDP_estimated(1:size(tau, 2));
        PDP_estimated = PDP_estimated / sum(PDP_estimated);
        FCF_estimated = fft([PDP_estimated; zeros(Links{BSID, UEID}.Channel.Nr.SamplesTotal - length(PDP_estimated), 1)]);
        FCF_estimated = circshift(FCF_estimated, [ceil(Links{BSID, UEID}.Channel.Nr.SamplesTotal / 2) 1]);

        % 计算 PDP 和 FCF 的均方误差（MSE）
        MSE_PDP = mean(abs(PDP_estimated - PDP_real).^2, 'all');
        MSE_FCF = mean(abs(FCF_estimated - FCF_real).^2, 'all');

        % 绘制估计和真实的 PDP
        figure(iBS*1000+iUE*100+10)
        set(gcf, 'Name', 'PDP对比')
        subplot(2,1,1)
        stem(tau, PDP_real, 'bx');
        title('Real Power Delay Profile (PDP)');
        xlabel('Delay');
        ylabel('Power');
        subplot(2,1,2)
        stem(tau, PDP_estimated, 'r--');
        title('Estimated Power Delay Profile (PDP)');
        xlabel('Delay');
        ylabel('Power');

        % 绘制估计和真实的 FCF
        figure(iBS*1000+iUE*100+11)
        set(gcf, 'Name', 'FCF对比')
        subplot(2,1,1)
        plot(frequency, abs(FCF_real));
        title('Real Frequency Correlation Function (FCF)');
        xlabel('Subcarrier');
        ylabel('Power');
        subplot(2,1,2)
        plot(frequency, abs(FCF_estimated), 'r--');
        title('Estimated Frequency Correlation Function (FCF)');
        xlabel('Subcarrier');
        ylabel('Power');

        % 显示 MSE
        disp(['MSE of PDP: ', num2str(MSE_PDP)]);
        disp(['MSE of FCF: ', num2str(MSE_FCF)]);

        % 统计带外泄露功率 OOBE
        fs = simParams.modulation.samplingRate;
        B = simParams.modulation.subcarrierSpacing * simParams.modulation.numerOfSubcarriers;
        N = length(Links{BSID, UEID}.ReceiveSignal);
        df = fs / N;
        psd = fft(Links{BSID, UEID}.ReceiveSignal);

        E = sum(abs(psd).^2);
        E_in = sum(abs(psd(1:ceil(B/df))).^2);
        E_oobe = E - E_in;
        disp(['OOBE: ', num2str(E_oobe/E)]);

        % 绘制发送信号的PSD三维图，横轴为OFDM符号，纵轴为子载波个数
        figure(iBS*1000+iUE*100+12)
        set(gcf, 'Name', '发送信号PSD');
        OFDMSymbol = floor(length(Links{BSID, UEID}.TransmitSignal) / simParams.modulation.numerOfSubcarriers);
        psd = fft(reshape(Links{BSID, UEID}.TransmitSignal(1:OFDMSymbol*simParams.modulation.numerOfSubcarriers), [], OFDMSymbol));
        mesh(10*log10(abs(psd)));
        title('Transmit Signal PSD');
        xlabel('OFDM Symbol');
        ylabel('Subcarrier');
        zlabel('Power/dB');

        % 绘制接收信号的PSD三维图，横轴为OFDM符号，纵轴为子载波个数
        figure(iBS*1000+iUE*100+13)
        set(gcf, 'Name', '接收信号PSD');
        OFDMSymbol = floor(length(Links{BSID, UEID}.UETotalSignal_) / simParams.modulation.numerOfSubcarriers);
        psd = fft(reshape(Links{BSID, UEID}.UETotalSignal_(1:OFDMSymbol*simParams.modulation.numerOfSubcarriers), [], OFDMSymbol));
        mesh(10*log10(abs(psd)));
        title('Receive Signal PSD');
        xlabel('OFDM Symbol');
        ylabel('Subcarrier');
        zlabel('Power/dB');
    end
end