%   Author
%       - Bashar Tahir, btahir@nt.tuwien.ac.at
%         (c) 2018 Institute of Telecommunications, TU Wien.
%         www.nt.tuwien.ac.at 

% Compile MEX C++ files in the "C++sourcefiles" folder. The compiled files
% are saved into the "+ChannelCodingComponents" folder.

fprintf('\n\nCompiling starts...\n');
fprintf('------------------------------------------------\n');
fprintf('Compiling Tail-Biting Convolutional decoder...\n');
mex '+Coding/C++sourcefiles/TBConvDecodeMEX.cpp' -output '+Coding/+ChannelCodingComponents/TBConvDecodeMEX';
fprintf('\nCompiling Turbo decoder...\n');
mex '+Coding/C++sourcefiles/turboDecodeMEX.cpp' -output '+Coding/+ChannelCodingComponents/turboDecodeMEX';
fprintf('\nCompiling LDPC decoder...\n');
mex '+Coding/C++sourcefiles/LDPCDecodeMEX.cpp' -output '+Coding/+ChannelCodingComponents/LDPCDecodeMEX';
fprintf('\nCompiling Polar decoder...\n');
mex '+Coding/C++sourcefiles/polarDecodeMEX.cpp' -output '+Coding/+ChannelCodingComponents/polarDecodeMEX';
fprintf('\nFinished compiling.\n');
fprintf('------------------------------------------------\n\n');