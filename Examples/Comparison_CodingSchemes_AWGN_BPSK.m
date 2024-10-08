%%
% BER Comparison Between Convolutional, Turbo, LDPC, and Polar Codes Using
% BPSK over the AWGN channel
%
% Author
%   Bashar Tahir, btahir@nt.tuwien.ac.at
%   (c) 2018 Institute of Telecommunications, TU Wien.
%   www.nt.tuwien.ac.at 

%% Initialization
clc; clear; close all;
cd('..');

% Assume that you need 1712 coded bits at the output of the encoder 
% (i.e. N = 1712)
CodeBlockLength = 1712;

% Set the code rate for all schemes
codeRate = 1/2;

ChannelCodingTBC = Coding.ChannelCoding(...
    'TB-Convolutional', ...          % Coding scheme
    'MAX-Log-MAP', ...               % Decoding algorithim
     codeRate ...                    % Code rate
);

ChannelCodingTurbo = Coding.ChannelCoding(...
    'Turbo', ...                     % Coding scheme
    'Linear-Log-MAP', ...               % Decoding algorithim
     codeRate, ...                   % Code rate
     8 ...                           % Decoding iterations
);

ChannelCodingLDPC = Coding.ChannelCoding(...
    'LDPC', ...                      % Coding scheme
    'PWL-Min-Sum', ...                   % Decoding algorithim
     codeRate, ...                   % Code rate
     16 ...                          % Decoding iterations
);

ChannelCodingPolar = Coding.ChannelCoding(...
    'Polar', ...                     % Coding scheme
    'CRC-List-SC', ...                   % Decoding algorithim
     codeRate, ...                   % Code rate
     8 ...                           % List size
);

% Use the "Update()" function to calculate how many bits should be input 
% to the encoder, and also prepare the object for that input

% Modulation Order
Qm = 1;

% Limited Soft Buffer: 1 for 100% buffer size (no limitation).
SoftBufferRatio = 1;

InputBitsLengthTBC = ChannelCodingTBC.update('Output', CodeBlockLength, codeRate, Qm, SoftBufferRatio);
InputBitsLengthTurbo = ChannelCodingTurbo.update('Output', CodeBlockLength, codeRate, Qm, SoftBufferRatio);
InputBitsLengthLDPC = ChannelCodingLDPC.update('Output', CodeBlockLength, codeRate, Qm, SoftBufferRatio);
InputBitsLengthPolar = ChannelCodingPolar.update('Output', CodeBlockLength, codeRate, Qm, SoftBufferRatio);

% This function must be called every time the required output length
% changes, or the code rate changes. Otherwise, it should be called only once
% outside the simulation repetition loop to reduce complexity.

NrRepetitions = 500;
SNRvecdB = linspace(0, 2, 5);

BER = zeros(5,length(SNRvecdB));

rng(6);

% Simulate the following coding schems
CodingSchemes = {'Uncoded', 'TBC', 'Turbo', 'LDPC', 'Polar'};
% To remove a coding scheme from the simulation, just remove it from
% the CodingSchemes cell.

%% Simulation
fprintf(['------- Started -------', '\n']);
st = tic;
for iiSNR = 1:length(SNRvecdB)
    stin = tic;
    SNR = SNRvecdB(iiSNR);
    for Codingi = 1:length(CodingSchemes)
        Coding = cell2mat(CodingSchemes(Codingi));
        if strcmp(Coding, 'Uncoded') == 0
            ChannelCoding1 = eval(['ChannelCoding' Coding]); 
            inputLength = eval(['InputBitsLength' Coding]);
        else
            ChannelCoding1 = 0;
            inputLength = ceil(CodeBlockLength*codeRate);
        end
        BERIn = 0;
        parfor Simulation = 1:NrRepetitions        
            % Generate the input bits
                InputBits = randi([0 1], inputLength, 1);
            % Encode the bits
            if strcmp(Coding, 'Uncoded') == 0
                CodedBits = ChannelCoding1.encode(InputBits); %#ok
            else
                CodedBits = InputBits;
            end
            % BPSK modulation
                y = CodedBits * -1;
                y(y==0) = 1;
            % AWGN channel
                sigma_n2 = 10^(-SNR/10);
                y = y + sqrt(sigma_n2) .* randn(length(y), 1);
            % LLRs calculation
                ChannelLLRs = 2.*y./sigma_n2;
            % Decode the bits
            if strcmp(Coding, 'Uncoded') == 0
                DecodedBits = ChannelCoding1.decode(ChannelLLRs);
            else
                DecodedBits = (1-sign(ChannelLLRs))./2;
            end 
            % BER mesurement
                BERIn = BERIn + sum(DecodedBits~=InputBits)/inputLength;   
        end
        BER(Codingi, iiSNR) = BERIn;
    end   
    fprintf(['Time remaining: ', num2str(round((toc(stin)*(length(SNRvecdB)-iiSNR)*100))/100), ' seconds', '\n']); 
end
toc(st)
BER = BER ./ NrRepetitions;

fprintf(['-------- Done --------', '\n']);

%% Plot results
semilogy(SNRvecdB, BER(1,:), 'linewidth', 1); hold on;
semilogy(SNRvecdB, BER(2,:), 'linewidth', 1); 
semilogy(SNRvecdB, BER(3,:), 'linewidth', 1); 
semilogy(SNRvecdB, BER(4,:), 'linewidth', 1); 
semilogy(SNRvecdB, BER(5,:), 'linewidth', 1); 

hold off;
xlabel('SNR [dB]');
ylabel('BER');
xlim([0 2.5]);
grid on;
legend('Uncoded', 'Convolutional (MAX-Log-MAP)','Turbo (Linear-Log-MAP, 8 iter)','LDPC (PWL-Min-Sum, 16 iter)','Polar (CRC-List-SC, 8 list)', 'location', 'northeast');
