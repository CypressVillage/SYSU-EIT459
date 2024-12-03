% 配置模型参数
cfgModel = winner2.wimparset;
cfgModel.NumTimeSamples = 1000; % 设置时间样本数
cfgModel.IntraClusterDsUsed = 'no'; % 不使用簇内延迟扩展
cfgModel.SampleDensity = 2e5; % 设置采样密度
cfgModel.PathLossModelUsed = 'yes'; % 使用路径损耗模型
cfgModel.ShadowingModelUsed = 'yes'; % 使用阴影模型

% 配置布局参数
BSAA = winner2.AntennaArray('UCA', 2 , 0.02); % 基站使用UCA-2阵列
MSAA = winner2.AntennaArray('ULA', 2, 0.01); % 移动站使用ULA-2阵列
cfgLayout = winner2.layoutparset(2, {1}, 1, [BSAA, MSAA], 200, 1);
cfgLayout.Stations(1).Pos(1:2) = [100; 100]; % 设置基站位置
cfgLayout.Stations(2).Pos(1:2) = [160; 100]; % 设置移动站位置
cfgLayout.PropagConditionVector = 1; % 设置传播条件为LOS

% 创建 WINNER II 信道对象
wimChan = comm.WINNER2Channel(cfgModel, cfgLayout);

% 生成信道系数
inputSignal = randn(1000, 2); % 输入信号
[outputSignal, channelCoefficients] = step(wimChan, inputSignal);

% 输出信道脉冲响应
disp(channelCoefficients);