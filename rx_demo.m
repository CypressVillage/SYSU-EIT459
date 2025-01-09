rng(42);  % 设置随机数种子为42
% %% USRP 参数初始化
radio = comm.SDRuReceiver(...
    'Platform', 'N200/N210/USRP2', ...
    'IPAddress', '192.168.10.2', ...
    'MasterClockRate', 100e6, ...
    'CenterFrequency', 1.5e9, ...
    'Gain', 10, ...
    'DecimationFactor', 10, ...
    'SamplesPerFrame', 5000 * 2.5, ... % 发送样点数 * 2.5
    'OutputDataType', 'double');

radio.OverrunOutputPort = true;

%% 相关输出参数初始化
len = uint32(0);
rcvdSignal = complex(zeros(20, 2));
% disp(simParams.usrpRX);
disp(['数据采集进行中...']);

hlog = dsp.SignalSink;

% 按照设定的接收次数循环接收
for idx = 1:20
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

% zc同步
transmitSignalLength = 4494; % 发送信号的长度
zcLength = 503;  % ZC序列的长度（可以根据需求调整）
zcSeed = 25;    % ZC序列的种子值（可以根据需求调整）
zcSequence = zadoffChuSeq(zcSeed, zcLength);  % 生成ZC序列

index = xcorr(zcSequence, UETotalSignal);
[~, maxIndex] = max(abs(index));
startID = length(UETotalSignal) - maxIndex + zcLength + 1;
endID = startID + transmitSignalLength - 1;
UETotalSignal = UETotalSignal(startID:endID);


save('RXusrp_data.mat', 'RXusrp_data', '-v7.3'); % 保存接收到的信息
