%% USRP 参数初始化
if strcmp(simParams.usrpRX.Platform, 'X310')
    radio = comm.SDRuReceiver(...
        'Platform', simParams.usrpRX.Platform, ...
        'IPAddress', simParams.usrpRX.Address, ...
        'IsTwinRXDaughterboard', simParams.usrpRX.IsTwinRXDaughterboard, ...
        'MasterClockRate', simParams.usrpRX.MasterClockRate, ...
        'ChannelMapping', simParams.usrpRX.ChannelMapping, ...
        'CenterFrequency', simParams.usrpRX.USRPCenterFrequency, ...
        'Gain', simParams.usrpRX.USRPGain, ...
        'DecimationFactor', simParams.usrpRX.USRPDecimationFactor, ...
        'SamplesPerFrame', simParams.usrpRX.USRPFrameLength, ...
        'OutputDataType', 'double');
else
    error(message('sdru:examples:UnsupportedPlatform', platform));
end

radio.OverrunOutputPort = true;

%% 相关输出参数初始化
len = uint32(0);
rcvdSignal = complex(zeros(simParams.usrpRX.USRPFrameLength, 2));
disp(simParams.usrpRX);
disp(['数据采集进行中...']);

hlog = dsp.SignalSink;

% 按照设定的接收次数循环接收
for idx = 1:simParams.usrpRX.numRxFrame
    % 当SDRu System对象输出有效时，继续访问SDRu系统对象输出，直到输出有效为止
    while len <= 0
        [rcvdSignal, len, overrun, timestamps] = step(radio);
    end
    if overrun ~= 0
        sprintf("接收数据未连续（不完整）\n");
    end
    % 当SDRu System对象输出有效时，解码接收到的消息
    if (len > 0)
        hlog(rcvdSignal);
    else
        sprintf("未接收到信号\n");
    end
    len = uint32(0);
end

disp(['数据采集结束！']);
release(radio);

RXusrp_data = hlog.Buffer;
save('RXusrp_data.mat', 'RXusrp_data', '-v7.3'); % 保存接收到的信息
