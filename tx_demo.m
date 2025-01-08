% USRP 初始化
radio = comm.SDRuTransmitter(...
    'Platform', 'N200/N210/USRP2', ...
    'IPAddress', '192.168.10.2', ...
    'MasterClockRate', 100e6, ...
    'CenterFrequency', 1.5e9, ...
    'Gain', 10, ...
    'InterpolationFactor', 20);

radio.UnderrunOutputPort = true;

tic; toc;

disp(radio);
% disp(simParams.usrpTX);
disp(['发送数据进行中...']);

txSig = load('TransmitSignal.mat', '-mat'); % 保存预先生成的发送信号txSig
txSig = txSig.var4_1;
% txSig = ones(10000, 1);

% 循环发送直至手动终止
while true
    radio(txSig);
    fprintf("dasdasda\n")
end

release(radio); % 释放USRP资源
