classdef UFMC < handle
    % =====================================================================    
    % Ljiljana Marijanovic, ljiljana.marijanovic@nt.tuwien.ac.at
    % (c) 2016 by Institute of Telecommunications, TU Wien
    % www.nt.tuwien.ac.at
    % =====================================================================    
    % This class represents Universal Filter Multicarrier (UFMC) system
    % model with corresponding transceiver parts. The modulation of data
    % symbols x and the demodulation of the received samples r are then 
    % performed by the methods ".Modulation(x)" and ".Demodulation(r)".
    % Also there are some useful functions such as
    % "PlotPowerSpectralDensity" and "PlotTransmitPower".
    % =====================================================================    
    
    properties (SetAccess = private)
        Nr   % for dimensionless parameters
        PHY  % for parameters with physical interpretation
        Implementation  % implmentation relevent parameters
    end
    
    methods
        function obj = UFMC(varargin)
            % Initialize parameters, set default values
            if numel(varargin)==10
                obj.Nr.Subcarriers = varargin{1};                           % number of subcarriers
                obj.Nr.MCSymbols = varargin{2};                             % number of multicarrier symbols
                obj.PHY.SubcarrierSpacing = varargin{3};                    % subcarrier spacing
                obj.PHY.SamplingRate = varargin{4};
                obj.PHY.IntermediateFrequency   = varargin{5};
                obj.PHY.TransmitRealSignal = varargin{6};
                obj.Nr.SubcarriersPerSubband = varargin{7};                 % number of subcarriers per one subband
                obj.PHY.ZeroPrefixLength = varargin{8};                     % length of zero prefix in samples
                obj.PHY.ZeroGuardTimeLength = varargin{9};                  % length of zero guard in seconds
                obj.PHY.FilterLength = varargin{10};                        % length of filter in samples (supposed to be larger for 1 sample than obj.PHY.ZeroPrefixLength)
            elseif numel(varargin)==0
                obj.Nr.Subcarriers = 12;
                obj.Nr.MCSymbols = 14;
                obj.PHY.SubcarrierSpacing = 15e3;
                obj.PHY.SamplingRate = obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing;
                obj.PHY.IntermediateFrequency = 0;
                obj.PHY.TransmitRealSignal = false;
                obj.Nr.SubcarriersPerSubband  = 12;
                obj.PHY.ZeroPrefixLength = 71;
                obj.PHY.ZeroGuardTimeLength = 0;
                obj.PHY.FilterLength = 1/(14*15e3);
            else
                error('Number of input variables must be either 0 (default values) or 10');
            end
            % calculate and set all dependent parameters
            obj.SetDependentParameters();
            
        end
        
        
        function SetDependentParameters(obj)
            % method that sets all parameters which are dependent on other
            % parameters
            
            if mod(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing),1)~=0
                obj.PHY.SubcarrierSpacing=obj.PHY.SamplingRate/(round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)));
                disp('Sampling Rate must be a multiple of the subcarrier spacing!');
                disp(['Therefore, the subcarrier spacing is set to: ' int2str(obj.PHY.SubcarrierSpacing) 'Hz']);
            end
            
            if (obj.PHY.SamplingRate<obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing)
                error('Sampling Rate must be higher: at least Number of Subcarriers times Subcarrier Spacing');
            end
            
            obj.Implementation.ZeroGuardSamples = round( obj.PHY.ZeroGuardTimeLength * obj.PHY.SamplingRate );
            obj.Implementation.FilterDuration = round(obj.PHY.FilterLength*obj.PHY.SamplingRate); % length of filter in samples
            obj.Implementation.ZeroPrefix = obj.PHY.ZeroPrefixLength/obj.PHY.SamplingRate; % length  of ZP in seconds
            obj.Implementation.IntermediateFrequency    = round( obj.PHY.IntermediateFrequency / obj.PHY.SubcarrierSpacing );
            obj.Implementation.FFTSize = obj.PHY.SamplingRate / obj.PHY.SubcarrierSpacing;
            obj.Implementation.TimeSpacing = round(obj.PHY.SamplingRate/obj.PHY.SubcarrierSpacing)+obj.PHY.ZeroPrefixLength;
            obj.Nr.Subbands =obj.Nr.Subcarriers /obj.Nr.SubcarriersPerSubband;

            if (rem(obj.Nr.Subbands,1)~=0)
                error('Number of subbands has to be integer; otherwise it is necessary to add zeros artificially');
            end
            obj.Nr.SamplesTotal = obj.Nr.MCSymbols *obj.Implementation.TimeSpacing + obj.Implementation.FilterDuration- 1;% total number of samples
            obj.PHY.dt = 1 / obj.PHY.SamplingRate;
            obj.PHY.TimeSpacing = obj.Implementation.TimeSpacing * obj.PHY.dt;
            obj.Implementation.NormalizationFactor = sqrt((obj.Implementation.FFTSize +obj.PHY.ZeroPrefixLength + obj.Implementation.FilterDuration))*sqrt(obj.Implementation.FFTSize/obj.Nr.Subcarriers);
            
        end
        
        
        % Set Functions
        function SetNrSubcarriers(varargin)
            % set the number of subcarriers
            
            obj = varargin{1};
            % set specific property
            obj.Nr.Subcarriers = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetNrMCSymbols(varargin)
            % set the number of symbols
            obj = varargin{1};
            % set specific property
            obj.Nr.MCSymbols = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetSubcarrierSpacing(varargin)
            % set the subcarrier spacing
            
            obj = varargin{1};
            % set specific property
            obj.PHY.SubcarrierSpacing = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetSamplingRate(varargin)
            % set the sampling rate
            
            obj = varargin{1};
            % set specific property
            obj.PHY.SamplingRate = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetTransmitRealSignal(varargin)
            % set real transmit signal indicator
            
            obj = varargin{1};
            % set specific property
            obj.PHY.TransmitRealSignal = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetSubbands(varargin)
            % set real transmit signal indicator
            
            obj = varargin{1};
            % set specific property
            obj.Nr.Subbands = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetSubcarriersPerRb(varargin)
            % set real transmit signal indicator
            
            obj = varargin{1};
            % set specific property
           obj.Nr.SubcarriersPerSubband = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetZeroPrefixLength(varargin)
            % set real transmit signal indicator
            
            obj = varargin{1};
            % set specific property
            obj.PHY.ZeroPrefixLength = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetFilterLength(varargin)
            % set real transmit signal indicator
            
            obj = varargin{1};
            % set specific property
            obj.Implementation.FilterDuration = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function txSignal = Modulation(varargin)
            % UFMC transmitter - Modulation of transmit data using
            % filtering. Total bandwidth is devided into the small subbands
            % and the modulation is carried out over each subband.
            % The input argument is a matrix
            % of size "Number of subcarriers" \times "Number of UFMC symbols (time)"
            % which represents the transmit data symbols.
            % The output is total transmit signal composed of
            % corresponding subband signals.
            
            obj = varargin{1};
            DataSymbols = varargin{2};
            
            DataSymbolsTemp = zeros(obj.Implementation.FFTSize,obj.Nr.MCSymbols);
            DataSymbolsTemp(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers),:) = DataSymbols.*obj.Implementation.NormalizationFactor;
            
            
            % hanning window
            f = hann(obj.Implementation.FilterDuration);
            
            if obj.Implementation.FilterDuration == 2
                % in this case we obtain two zero elements that causes
                % problem with normalization
                f = [1;0];
            end
            
            % normalization
            f = f./norm(f);
            
            totalTxSignal = 0;
            for m = 1:obj.Nr.Subbands
                % subband data
                n = ((m-1)* obj.Nr.SubcarriersPerSubband +1:(m-1)*obj.Nr.SubcarriersPerSubband +obj.Nr.SubcarriersPerSubband)+obj.Implementation.IntermediateFrequency ;
                b = reshape(DataSymbolsTemp(n,:),obj.Nr.SubcarriersPerSubband ,obj.Nr.MCSymbols);
                
                % center frequency of each subband
                centralFreq = ((obj.Nr.SubcarriersPerSubband +1)/2 + (m-1)*obj.Nr.SubcarriersPerSubband) + obj.Implementation.IntermediateFrequency;
                
                for l = 1:obj.Implementation.FilterDuration
                    freqShift (l) = exp(2*pi*1i*(l-1)*(centralFreq-1)/obj.Implementation.FFTSize);
                end
                
                % filter in time domain
                fTime = f'.*freqShift;
                % filter in frequency domain
                fFreq = fft(fTime.',obj.Implementation.FFTSize);
                
                % pre-equalization part
                if obj.Implementation.FilterDuration == 1
                    fFreq = fFreq';
                end
                b1 = b./repmat(fFreq(n),1,obj.Nr.MCSymbols);
                b1 = [zeros(obj.Implementation.IntermediateFrequency + obj.Nr.SubcarriersPerSubband *(m-1),obj.Nr.MCSymbols);b1;zeros(obj.Nr.SubcarriersPerSubband*(obj.Nr.Subbands-m),obj.Nr.MCSymbols)];
                
                % IFFT part (time-domain signal)
                txTime = ifft(b1,obj.Implementation.FFTSize);
                
                % Zero-padding part
                txTimeZP = [txTime;zeros(obj.PHY.ZeroPrefixLength,obj.Nr.MCSymbols)];
                
                % subband-filtering part
                
                txFiltered = conv(txTimeZP(:),fTime);

                % addition of signals from multiple subbands
                totalTxSignal = totalTxSignal + txFiltered;
            end
            txSignal = totalTxSignal;
        end
        
        
        function rxSignal = Demodulation(varargin)
            % UFMC receiver - Demodulation of receive signal
            % The input argument is a signal after the channel 
            % currupted also by the noise.
            % The output is a signal after FFT. (frequency domain)
            % If there is a receive filtering, "rxSignal" is one obtained
            % after frequency domain processing (i.e., equalization)
            
            obj = varargin{1};
            ReceivedSignal = varargin{2};
            % In order to appropriately reshape our input signal taking
            % into account exact number of the UFMC symbols this step must
            % be done
            ReceivedSignal = ReceivedSignal(1:end-obj.Implementation.FilterDuration+1);
            ReceivedSignalResh = reshape(ReceivedSignal,[] ,obj.Nr.MCSymbols);
            
            % creating the tail, see R1-165014
            guardLength = size(ReceivedSignalResh,1)-obj.Implementation.FFTSize;
            receiveSignalTemp = ReceivedSignalResh(1:obj.Implementation.FFTSize,:);
            
            % Copying of tail at the beginning of signal, see R1-165014
            receiveSignalTemp(1:guardLength,:) =  ReceivedSignalResh(1:guardLength,:) +  ReceivedSignalResh(obj.Implementation.FFTSize+1:end,:);
            
            % N-FFT part
            receiveSignalFFT  = fft(receiveSignalTemp,obj.Implementation.FFTSize);
            rxSignal = receiveSignalFFT(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers),:)/obj.Implementation.NormalizationFactor;
        end
        
        function Pn = GetSymbolNoisePower(varargin)
            % calculates the symbol noise power (after demodulation)
            % The input argument is the noise power in the time domain. 
            obj = varargin{1};
            Pn_time = varargin{2};
            Pn = Pn_time .*(obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing/(obj.PHY.SamplingRate));
        end
        
        function TimePos = GetTimeIndexMidPos(obj)
            % returns a vector which represents the discete time position
            % (middle positions)
            TimePos = obj.Implementation.ZeroGuardSamples+obj.PHY.ZeroPrefixLength+obj.Implementation.FilterDuration+round(obj.Implementation.FFTSize/2)+1+ (0:obj.Nr.MCSymbols-1)*obj.Implementation.TimeSpacing;
        end
        
        function [TransmitPower,Time] = PlotTransmitPower(obj, Rx)
            % plot the expected transmit power over time. The input 
            % argument represents the correlation of the data symbols. 
            % If no input argument is specified, an identity matrix is 
            % assumed (uncorrelated data) 
            
            if exist('Rx','var')
                [V,D] = eig(Rx);
            else 
                % assume that Rx is an identity matrix, that is,
                % uncorrelated symbols
                V = eye(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
                D = V;               
            end
            D=sqrt(D);
            TransmitPower = zeros(obj.Nr.SamplesTotal,1);
            for i_lk = 1:obj.Nr.Subcarriers*obj.Nr.MCSymbols
                % For OFDM.PHY.TransmitRealValuedSignal == 0, the modulation is not linear => we have to consider 1 and +j.
                TransmitPower = TransmitPower+(abs(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2+abs(obj.Modulation(reshape(1j*V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2)/2;
                if mod(i_lk,1000)==0
                    disp([int2str(i_lk/(obj.Nr.Subcarriers*obj.Nr.MCSymbols )*100) '%']);
                end
            end
            Time = (0:length(TransmitPower)-1)*obj.PHY.dt;
            if nargout==0
                plot(Time,TransmitPower);
                ylabel('Transmit Power');
                xlabel('Time(s)');
            end
        end
        
        function [PowerSpectralDensity,Frequency] = PlotPowerSpectralDensity(obj,Rx)
            % plot the power spectral density. The input argument 
            % represents the correlation of the data symbols. If no input
            % argument is specified, an identity matrix is assumed 
            % (uncorrelated data) 
            
            if exist('Rx','var')
                [V,D] = eig(Rx);
            else 
                V = eye(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
                D = V;               
            end
            D=sqrt(D);
            PowerSpectralDensity = zeros(obj.Nr.SamplesTotal,1);
            for i_lk = 1:obj.Nr.Subcarriers*obj.Nr.MCSymbols
                PowerSpectralDensity = PowerSpectralDensity+abs(fft(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk))).^2;
                if mod(i_lk,1000)==0
                    disp([int2str(i_lk/(obj.Nr.Subcarriers*obj.Nr.MCSymbols )*100) '%']);
                end
            end
            Frequency = (0:length(PowerSpectralDensity)-1)*1/(length(PowerSpectralDensity)*obj.PHY.dt);
            PowerSpectralDensity=PowerSpectralDensity/length(PowerSpectralDensity)^2/Frequency(2)^2;
            if nargout==0
                plot(Frequency,10*log10(PowerSpectralDensity));
                ylabel('Power Spectral Density (dB/Hz)');
                xlabel('Frequency (Hz)');
            end
        end
        
        
    end
end



