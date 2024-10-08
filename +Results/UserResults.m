classdef UserResults < handle
    % this class contains all user specific simulation results
    
    properties
        nFrames                 % number of simulated frames
        nSweep                  % number of simulated sweep values
        SNR                     % Signal to Noise Ratio (per frame and sweep value)
        nCodeBits               % number of code bits (per frame and sweep value, for each spatial stream)
        nDataBits               % number of data bits (per frame and sweep value, for each spatial stream)
        nSpatialStreamsUsed     % number of employed spatial streams (per frame and sweep value)
        mcsUsed                 % employed Modulation and Coding Scheme (per frame and sweep value, for each spatial stream)
        bitErrorsCoded          % number of bit errors calculated from the data bits (per frame and sweep value, for each spatial stream)
        bitErrorsUncoded        % number of bit errors calculated from the code bits (per frame and sweep value, for each spatial stream)
        frameErrors             % frame error per frame and sweep value (per frame and sweep value, for each spatial stream)
        BERCoded                % Bit Error Ratio calculated from coded bit errors (per sweep value, for each spatial stream, with confidence)
        BERUncoded              % Bit Error Ratio calculated from uncoded bit errors (per sweep value, for each spatial stream, with confidence)
        FER                     % Frame Error Ratio (per sweep value, with confidence)
        throughput              % throughput calculated from error free frames, for all spatial streams combined (per sweep value, with confidence)
        channelMSE              % channel estimation mean squared error
        PAPR                    % Peak-to-Average Power Ratio
    end
    
    methods
        function obj = UserResults( nFrames, nSweep )
            % user results class constructor
            
            if nargin > 0
                obj.nFrames = nFrames;
                obj.nSweep  = nSweep;

                % initialize result parameters
                obj.SNR                                 = NaN( nFrames, nSweep );
                obj.nSpatialStreamsUsed.values          = NaN( nFrames, nSweep );
                obj.mcsUsed.values                      = NaN( nFrames, nSweep );
                obj.PAPR.values                         = NaN( nFrames, nSweep);
                obj.PAPR.CDF                            = NaN( nFrames*nSweep + 1, 1);
                obj.PAPR.DataPoints                     = NaN( nFrames*nSweep + 1, 1);
                
                % initialize results that have per stream values
                obj.nCodeBits.values                    = NaN( nFrames, nSweep );   % total
                obj.nCodeBits.valuesPerStream           = cell( nFrames, nSweep );  % per stream
                obj.nDataBits.values                    = NaN( nFrames, nSweep );   % total
                obj.nDataBits.valuesPerStream           = cell( nFrames, nSweep );  % per stream
                obj.bitErrorsCoded.values               = NaN( nFrames, nSweep );   % total
                obj.bitErrorsCoded.valuesPerStream      = cell( nFrames, nSweep );  % per stream
                obj.bitErrorsUncoded.values             = NaN( nFrames, nSweep );   % total
                obj.bitErrorsUncoded.valuesPerStream    = cell( nFrames, nSweep );  % per stream
                obj.frameErrors.values                  = NaN( nFrames, nSweep );   % total
                obj.frameErrors.valuesPerStream         = cell( nFrames, nSweep );  % per stream        

                % initialize result parameters that have confidence
                obj.throughput.values           = NaN( nFrames, nSweep );
                obj.throughput.mean             = NaN( 1, nSweep );
                obj.throughput.confidence       = NaN( 2, nSweep );
                
                obj.BERCoded.values             = NaN( nFrames, nSweep );
                obj.BERCoded.mean               = NaN( 1, nSweep );
                obj.BERCoded.confidence         = NaN( 2, nSweep );
                
                obj.BERUncoded.values           = NaN( nFrames, nSweep );
                obj.BERUncoded.mean             = NaN( 1, nSweep );
                obj.BERUncoded.confidence       = NaN( 2, nSweep );
                
                obj.FER.values                  = NaN( nFrames, nSweep );
                obj.FER.mean                    = NaN( 1, nSweep );
                obj.FER.confidence              = NaN( 2, nSweep );
                
                obj.channelMSE.values           = NaN( nFrames, nSweep );
                obj.channelMSE.mean             = NaN( 1, nSweep );
                obj.channelMSE.confidence       = NaN( 2, nSweep );
            end
            
        end
    end
    
end

