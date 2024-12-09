%% 调用WINNER-II信道模型
ScenarioVectors = [1 2 3 5 10]; %1=A1, 2=A2, 3=B1, 4=B2, 5=B3, 6=B4, 10=C1, 11=C2, 12=C3, 13=C4, 14=D1, 15=D2a
PropagConditionVectors = [1 0 1 1 0]; %LOS = 1 and NLOS = 0

cfgWim = winner2.wimparset;
samplingRate = 15.36e6;
cfgWim.CenterFrequency = 2.5e9; % 设置中心频率
cfgWim.DelaySamplingInterval = 1/samplingRate;
BSAA = winner2.AntennaArray();
MSAA = winner2.AntennaArray();
MSIdx = [1];
BSIdx = {1};
NL = 1; % 单径
rndSeed = 5;
layoutpar = winner2.layoutparset(MSIdx,BSIdx,NL,[BSAA,MSAA],[],rndSeed);

for iChannel = 1:length(ScenarioVectors)
    layoutpar.ScenarioVector = ScenarioVectors(iChannel);
    layoutpar.PropagConditionVector = PropagConditionVectors(iChannel);
    [H,pathdealy,~] = winner2.wim(cfgWim,layoutpar);
    [~, ~, row, col] = size(H{1,1});
    A = zeros(row, col);
    % 取出信道矩阵
    for ii = 1:row
        for jj = 1:col
            A(ii, jj) = H{1,1}(1,1,ii,jj);
        end
    end
    PDP = mean(A, 2);
    PDP = abs(PDP).^2;
    PDP_mat = zeros(size(PDP));
    % 将同一条时延路径上的能量相加
    for index = 1:length(pathdealy)
        PDP_mat(index) = sum(PDP(pathdealy==pathdealy(index)));
    end

    figure(2000 + iChannel)
    stem(pathdealy, PDP_mat)

    xlabel('delay[s]')
    ylabel('Power Delay Profile')
end