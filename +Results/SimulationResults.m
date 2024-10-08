classdef SimulationResults < handle
    % this class contains and processes the simulation results
    
    properties
        userResults     % per user simulation results
        nFrames         % number of simulated frames
        nSweep          % number of sweep values
        nUE             % number of users per base station
        nBS             % number of base stations
        frameDuration   % time duration of a frame
        resType         % results type {downlink, uplink, D2D}
    end

    methods
        function obj = SimulationResults( nFrames, nSweep, nBS, nUE, time, resType )
            % simulationResults class constructor
            
            if nargin > 0
                obj.nFrames         = nFrames;
                obj.nSweep          = nSweep;
                obj.nUE             = nUE;
                obj.nBS             = nBS;
                obj.frameDuration   = time;
                
                % initialize results type
                switch resType
                    case {'downlink', 'uplink', 'D2D'}
                        obj.resType = resType;
                    otherwise
                        error('Results type unkown.');
                end
                
                % initialize user results
                obj.userResults = Results.UserResults( nFrames, nSweep );
                for iUE = 1:nUE
                    obj.userResults(iUE) = Results.UserResults( nFrames, nSweep );
                end
            end
            
        end
        
        function collectResults( obj, sweepResults, UE )
            % this method collects temporary results from parallel
            % procesing
            
            switch obj.resType
                case {'downlink'}
                    for iSweep = 1: obj.nSweep
                        tmpResults = sweepResults{iSweep};
                        for iUE = 1:obj.nUE
                            for iFrame = 1:obj.nFrames
                                obj.userResults(iUE).SNR(iFrame,iSweep)                                 = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.SNR;
                                obj.userResults(iUE).nSpatialStreamsUsed.values(iFrame,iSweep)          = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.nStreams;
                                
                                % results per stream and for all streams in total
                                obj.userResults(iUE).nDataBits.valuesPerStream{iFrame,iSweep}           = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.nDataBits;
                                obj.userResults(iUE).nDataBits.values(iFrame,iSweep)                    = sum( obj.userResults(iUE).nDataBits.valuesPerStream{iFrame,iSweep} );
                                obj.userResults(iUE).nCodeBits.valuesPerStream{iFrame,iSweep}           = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.nCodeBits;
                                obj.userResults(iUE).nCodeBits.values(iFrame,iSweep)                    = sum( obj.userResults(iUE).nCodeBits.valuesPerStream{iFrame,iSweep} );

                                obj.userResults(iUE).bitErrorsUncoded.valuesPerStream{iFrame,iSweep}    = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.nBitError;
                                obj.userResults(iUE).bitErrorsUncoded.values(iFrame,iSweep)             = sum( obj.userResults(iUE).bitErrorsUncoded.valuesPerStream{iFrame,iSweep} );
                                obj.userResults(iUE).bitErrorsCoded.valuesPerStream{iFrame,iSweep}      = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.nCodedBitError;
                                obj.userResults(iUE).bitErrorsCoded.values(iFrame,iSweep)               = sum( obj.userResults(iUE).bitErrorsCoded.valuesPerStream{iFrame,iSweep} );
                                
                                obj.userResults(iUE).frameErrors.valuesPerStream{iFrame,iSweep}         = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.frameError;
                                obj.userResults(iUE).frameErrors.values(iFrame,iSweep)                  = mean( obj.userResults(iUE).frameErrors.valuesPerStream{iFrame,iSweep} );
                                
                                obj.userResults(iUE).channelMSE.values(iFrame,iSweep)                   = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.channelMSE;
                                
                                obj.userResults(iUE).PAPR.values(iFrame,iSweep)                         = tmpResults{UE{iUE}.TransmitBS(1), UE{iUE}.ID, iFrame}.PAPR;
                            end %for iFrame
                        end %for iUE
                    end %for iSweep
                case {'uplink'}
                    for iSweep = 1: obj.nSweep
                        tmpResults = sweepResults{iSweep};
                        for iUE = 1:obj.nUE
                            for iFrame = 1:obj.nFrames
                                obj.userResults(iUE).SNR(iFrame,iSweep)                                 = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.SNR;
                                obj.userResults(iUE).nSpatialStreamsUsed.values(iFrame,iSweep)          = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.nStreams;
                                
                                % results per stream and for all streams in total
                                obj.userResults(iUE).nDataBits.valuesPerStream{iFrame,iSweep}           = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.nDataBits;
                                obj.userResults(iUE).nDataBits.values(iFrame,iSweep)                    = sum( obj.userResults(iUE).nDataBits.valuesPerStream{iFrame,iSweep} );
                                obj.userResults(iUE).nCodeBits.valuesPerStream{iFrame,iSweep}           = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.nCodeBits;
                                obj.userResults(iUE).nCodeBits.values(iFrame,iSweep)                    = sum( obj.userResults(iUE).nCodeBits.valuesPerStream{iFrame,iSweep} );
                                
                                obj.userResults(iUE).bitErrorsUncoded.valuesPerStream{iFrame,iSweep}    = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.nBitError;
                                obj.userResults(iUE).bitErrorsUncoded.values(iFrame,iSweep)             = sum( obj.userResults(iUE).bitErrorsUncoded.valuesPerStream{iFrame,iSweep} );
                                obj.userResults(iUE).bitErrorsCoded.valuesPerStream{iFrame,iSweep}      = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.nCodedBitError;
                                obj.userResults(iUE).bitErrorsCoded.values(iFrame,iSweep)               = sum( obj.userResults(iUE).bitErrorsCoded.valuesPerStream{iFrame,iSweep} );

                                obj.userResults(iUE).frameErrors.valuesPerStream{iFrame,iSweep}         = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.frameError;
                                obj.userResults(iUE).frameErrors.values(iFrame,iSweep)                  = mean( obj.userResults(iUE).frameErrors.valuesPerStream{iFrame,iSweep} );
                                
                                obj.userResults(iUE).channelMSE.values(iFrame,iSweep)                   = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.channelMSE;
                                
                                obj.userResults(iUE).PAPR.values(iFrame,iSweep)                         = tmpResults{UE{iUE}.ID, UE{iUE}.ReceiveBS(1), iFrame}.PAPR;
                            end %for iFrame
                        end %for iUE
                    end %for iSweep
                otherwise
                    error('Results type unkown.');
            end %switch obj.resType
            
        end %function collectResults
        
        function postProcessResults( obj )
            % this function post processes the saved simulation results,
            % i.e., it calculates mean values and confidence intervals
            
            % function handles for bootci
            ber_func    = @(bit_error,nr_of_bits) sum(bit_error)./sum(nr_of_bits);
            mean_func   = @(inputValue) mean(inputValue);
            
            % calculate confidence intervals
            for iUE = 1:obj.nUE
                % coded Bit Error Ratio
                obj.userResults(iUE).BERCoded.values          = obj.userResults(iUE).bitErrorsCoded.values ./ obj.userResults(iUE).nDataBits.values;
                obj.userResults(iUE).BERCoded.mean            = nanmean( obj.userResults(iUE).bitErrorsCoded.values ./ obj.userResults(iUE).nDataBits.values );
                mask = ~isnan( obj.userResults(iUE).bitErrorsCoded.values );
                if mask
                    bitErrosCodedMasked = reshape( obj.userResults(iUE).bitErrorsCoded.values( mask ), [], obj.nSweep);
                    nDataBitsMasked    = reshape( obj.userResults(iUE).nDataBits.values, [], obj.nSweep);
                    obj.userResults(iUE).BERCoded.confidence  = bootci( 2000, {ber_func, bitErrosCodedMasked, nDataBitsMasked }, 'alpha', 0.05);
                    obj.userResults(iUE).BERCoded.confidence  = repmat([-1;1], 1, obj.nSweep) .* ( obj.userResults(iUE).BERCoded.confidence - repmat(obj.userResults(iUE).BERCoded.mean, 2, 1) );
                else
                    obj.userResults(iUE).BERCoded.confidence  = NaN(2, obj.nSweep );
                end

                % uncoded Bit Error Ratio
                obj.userResults(iUE).BERUncoded.values          = obj.userResults(iUE).bitErrorsUncoded.values ./ obj.userResults(iUE).nCodeBits.values;
                obj.userResults(iUE).BERUncoded.mean            = nanmean( obj.userResults(iUE).bitErrorsUncoded.values ./ obj.userResults(iUE).nCodeBits.values );
                mask = ~isnan( obj.userResults(iUE).bitErrorsUncoded.values );
                if mask
                    bitErrosUncodedMasked   = reshape( obj.userResults(iUE).bitErrorsUncoded.values( mask ), [], obj.nSweep);
                    nBitsCodedMasked      = reshape( obj.userResults(iUE).nCodeBits.values, [], obj.nSweep);
                    obj.userResults(iUE).BERUncoded.confidence  = bootci( 2000, {ber_func, bitErrosUncodedMasked, nBitsCodedMasked }, 'alpha', 0.05);
                    obj.userResults(iUE).BERUncoded.confidence  = repmat([-1;1], 1, obj.nSweep) .* ( obj.userResults(iUE).BERUncoded.confidence - repmat(obj.userResults(iUE).BERUncoded.mean, 2, 1) );
                else
                    obj.userResults(iUE).BERUncoded.confidence  = NaN(2, obj.nSweep );
                end
 
                % Frame Error Ratio
                framesSent = double( ~isnan( obj.userResults(iUE).nDataBits.values ) );
                obj.userResults(iUE).FER.values = obj.userResults(iUE).frameErrors.values ./ framesSent;
                obj.userResults(iUE).FER.mean = mean( obj.userResults(iUE).frameErrors.values ./ framesSent );
                mask = ~isnan( obj.userResults(iUE).FER.values );
                if mask
                    frameErrorsMasked  = reshape( obj.userResults(iUE).frameErrors.values( mask ), [], obj.nSweep);
                    nFramesMasked      = reshape( framesSent( mask ), [], obj.nSweep);
                    obj.userResults(iUE).FER.confidence = bootci( 2000, {ber_func, frameErrorsMasked, nFramesMasked }, 'alpha', 0.05);
                    obj.userResults(iUE).FER.confidence = repmat([-1;1], 1, obj.nSweep) .* ( obj.userResults(iUE).FER.confidence - repmat(obj.userResults(iUE).FER.mean, 2, 1) );
                else
                    obj.userResults(iUE).FER.confidence = NaN(2, obj.nSweep );
                end
                
                % throughput
                obj.userResults(iUE).throughput.values = ( obj.userResults(iUE).nDataBits.values .* (1-double(obj.userResults(iUE).frameErrors.values)) ) / obj.frameDuration;
                obj.userResults(iUE).throughput.mean = mean( obj.userResults(iUE).throughput.values );
                mask = ~isnan( obj.userResults(iUE).throughput.values );
                if mask
                    throughputMasked  = reshape( obj.userResults(iUE).throughput.values( mask ), [], obj.nSweep);
                    obj.userResults(iUE).throughput.confidence = bootci( 2000, {mean_func, throughputMasked }, 'alpha', 0.05);
                    obj.userResults(iUE).throughput.confidence = repmat([-1;1], 1, obj.nSweep) .* ( obj.userResults(iUE).throughput.confidence - repmat(obj.userResults(iUE).throughput.mean, 2, 1) );
                else
                    obj.userResults(iUE).throughput.confidence = NaN(2, obj.nSweep );
                end
                
                % channel estimation error
                obj.userResults(iUE).channelMSE.mean = mean(obj.userResults(iUE).channelMSE.values);
                mask = ~isnan( obj.userResults(iUE).channelMSE.values );
                if mask
                    channelMSEMasked = reshape( obj.userResults(iUE).channelMSE.values( mask ), [], obj.nSweep);
                    obj.userResults(iUE).channelMSE.confidence = bootci( 2000, {mean_func, channelMSEMasked }, 'alpha', 0.05);
                    obj.userResults(iUE).channelMSE.confidence = repmat([-1;1], 1, obj.nSweep) .* ( obj.userResults(iUE).channelMSE.confidence - repmat(obj.userResults(iUE).channelMSE.mean, 2, 1) );
                else
                    obj.userResults(iUE).channelMSE.confidence = NaN(2, obj.nSweep );
                end
                
                % PAPR CDF calculation
                [obj.userResults(iUE).PAPR.CDF, obj.userResults(iUE).PAPR.DataPoints] = ecdf(obj.userResults(iUE).PAPR.values(:));

            end %for iUE

        end %function postProcessResults
             
    end
    
end

