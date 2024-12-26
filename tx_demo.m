% USRP 初始化
if strcmp(simParams.usrpTX.Platform, 'N200')
    radio = comm.SDRuTransmitter(...
        'Platform', simParams.usrpTX.Platform, ...
        'IPAddress', simParams.usrpTX.Address, ...
        'MasterClockRate', simParams.usrpTX.MasterClockRate, ...
        'ChannelMapping', simParams.usrpTX.ChannelMapping, ...
        'CenterFrequency', simParams.usrpTX.USRPCenterFrequency, ...
        'Gain', simParams.usrpTX.USRPGain, ...
        'InterpolationFactor', simParams.usrpTX.USRPInterpolationFactor);
else
    error(message('sdru:examples:UnsupportedPlatform', platform));
end

radio.UnderrunOutputPort = true;

tic; toc;

disp(radio);
disp(simParams.usrpTX);
disp(['发送数据进行中...']);

save('txSig.mat', 'txSig'); % 保存预先生成的发送信号txSig

% 循环发送直至手动终止
while true
    radio(txSig);
end

release(radio); % 释放USRP资源
