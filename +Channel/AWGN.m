function noise = AWGN( noisePower, totalSamples, nAntennas )
% Generate complex Additive White Gaussian Noise

    noise = sqrt(noisePower/2) * complex(randn(totalSamples, nAntennas), randn(totalSamples, nAntennas));
end

