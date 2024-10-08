%%
% Short block length comparison between Convolutional, Turbo, LDPC, and Polar
% Codes.
%
% Author
%   Bashar Tahir, btahir@nt.tuwien.ac.at
%   (c) 2018 Institute of Telecommunications, TU Wien.
%   www.nt.tuwien.ac.at 
%%
clc;clear;
cd('..');

%% Simulation setup
SNRvecdB = -7:0.5:-1;  % SNR range
numR = 1000;           % number of simulated frames

%% Information length and code rate
K = 48;
R = 1/6;
R_LDPC = 0.256; % For LDPC, we use a higher code rate in order to compensate for the rate reduction caused by the addition of filler bits.

%% Modulation
modOrder = 4;
Mapper = Modulation.SignalConstellation(modOrder, 'QAM', 0); 

%% The channel coders
coderLDPC = Coding.ChannelCoding('LDPC', 'PWL-Min-Sum', R_LDPC, 32);  
coderPolar = Coding.ChannelCoding('Polar', 'CRC-List-SC', R, 32);
coderTurbo = Coding.ChannelCoding('Turbo', 'Linear-Log-MAP', R, 16);      
coderConv = Coding.ChannelCoding('TB-Convolutional', 'MAX-Log-MAP', R);      

%% Initialization
% Force CRC-16
coderLDPC.CRCPoly = '16';
coderLDPC.CRCLength = 16;
coderPolar.CRCPoly = '16';
coderPolar.CRCLength = 16;
coderTurbo.CRCPoly = '16';
coderTurbo.CRCLength = 16;
coderConv.CRCPoly = '16';
coderConv.CRCLength = 16;

% Update the objects
coderLDPC.update('Input', K, R_LDPC, log2(modOrder), 1);
coderPolar.update('Input', K, R, log2(modOrder), 1);
coderTurbo.update('Input', K, R, log2(modOrder), 1);
coderConv.update('Input', K, R, log2(modOrder), 1);

%% Simulation loop
Results = zeros(numR, length(SNRvecdB), 4);
fprintf(['------- Started -------', '\n']);
d = tic;
parfor m = 1:numR
    st = tic;
    ind = 1;
    resultsPerCoder = zeros(length(SNRvecdB), 4);
    for coder = [coderConv, coderTurbo, coderLDPC, coderPolar]
        Resultsin = zeros(1, length(SNRvecdB));
        for iSNR = 1:length(SNRvecdB)
            SNRlinear = 10^(SNRvecdB(iSNR)/10);   
            sigma_n2 = 1/SNRlinear;
            
            InputBits = randi([0 1], K, 1);
            CodedBits = coder.encode(InputBits);
            padded = zeros(ceil(length(CodedBits)/2)*2-length(CodedBits),1);
            yCoded = Mapper.Bit2Symbol([CodedBits;padded]);

            yCoded = yCoded + sqrt(sigma_n2/2) .* complex(randn(length(yCoded),1), randn(length(yCoded),1));

            ChannelLLRsCoded = -Mapper.LLR_AWGN(yCoded, sigma_n2);
            ChannelLLRsCoded(end-length(padded)+1:end) = [];
            DecodedBits = coder.decode(ChannelLLRsCoded);

%             Resultsin(iSNR) =  sum(DecodedBits~=InputBits)/(K);          % uncomment this to get BER results.
            Resultsin(iSNR) = double(sum(DecodedBits~=InputBits)>0);       % uncomment this to get ideal FER results.
%             Resultsin(iSNR) = Coder.CRCDetectionResult;                  % uncomment this to get FER results based on the CRC calculation.
        end
        resultsPerCoder(:,ind) = Resultsin;
        ind = ind + 1;
    end  
    Results(m,:,:) = resultsPerCoder;
    en = toc(st);
    fprintf('Repetition %i/%i.\n',m,numR); 
end
fprintf(['-------- Done --------', '\n']);
toc(d)
%% Plot results
% uncoded AWGN
semilogy(SNRvecdB, 0.5*erfc(sqrt(0.5*10.^(SNRvecdB./10))), 'r', 'linewidth', 1);  hold on;

% coded
ind = 1;
for coder = [coderConv, coderTurbo, coderLDPC, coderPolar]
    proResults = Results(:,:,ind);
    ind = ind + 1;
    fprintf(['Calculating confidence intervals...', '\n']);
    codewordsSent = double(ones(size(proResults)));
    values = proResults ./ codewordsSent;
    meanVal = mean(values ./ codewordsSent);
    mask = ~isnan(values);
    if mask
        ErrorsMasked  = reshape(values( mask ), [], length(SNRvecdB));
        nMasked      = reshape( codewordsSent( mask ), [], length(SNRvecdB));
        conf = bootci( 2000, {@(error,nrError) sum(error)./sum(nrError), ErrorsMasked, nMasked }, 'alpha', 0.05);
        conf = repmat([-1;1], 1, length(SNRvecdB)) .* ( conf - repmat(meanVal, 2, 1) );
    else
        conf = NaN(2, length(SNRvec));
    end
    fprintf(['Done.', '\n']);
    errorbar(SNRvecdB, meanVal, conf(1,:),conf(2,:));
end
% plot settings
hold off;
xlabel('SNR (dB)');
ylabel('FER');
grid on;
ylim([1e-4 1e0]);
legend('uncoded', 'convolutional', 'turbo', 'LDPC', 'polar','location','best')
cd('Examples')