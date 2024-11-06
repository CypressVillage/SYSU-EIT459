classdef FastFading < handle
    % =====================================================================        
    % This MATLAB class represents an implementation of a time-variant
    % channel. The channel parameters are initialized by the class
    % constructor. The convolution is performed by the method
    % ".Convolution(.)". A new channel realization is obtained by 
    % ".NewRealization;".
    % Additionally, there exist some other useful methods, such as
    % ".GetRmsDelaySpread" or ".PlotPowerDelayProfile"
    % =====================================================================    
    % Ronald Nissel, rnissel@nt.tuwien.ac.at
    % (c) 2017 by Institute of Telecommunications, TU Wien
    % www.tc.tuwien.ac.at
    % =====================================================================    
    
    
    properties (SetAccess = private)
        PHY
        Nr
        Implementation
        ImpulseResponse
    end
    
    methods
        function obj = FastFading(...
                SamplingRate,...                                    % Sampling rate (Samples/s)
                PowerDelayProfile,...                               % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
                SamplesTotal,...                                    % Number of total samples
                MaximumDopplerShift,...                             % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
                DopplerModel,...                                    % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
                Paths,...                                           % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum
                timeCorrelation,...                                 % whether the channel is correlated in time (across frames) or not
                spatialCorrelation,...                              % spatial correlation of MIMO channels
                correlationCoeffTx,...                              % tansmit side correlation coefficient
                correlationCoeffRx,...                              % receive side correlation coefficient
                nTxAntennas,...                                     % Number of transmit antennas
                nRxAntennas,...                                     % Number of receive antennas
                WarningIfSampleRateDoesNotMatch,...                 % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
                K,...                                               % TWDP model parameter K
                Delta ...                                            % TWDP model parameter delta
                )
            
            % Initialize parameters
            obj.PHY.SamplingRate                    = SamplingRate;
            obj.Nr.SamplesTotal                     = SamplesTotal;
            obj.Nr.txAntennas                       = nTxAntennas;
            obj.Nr.rxAntennas                       = nRxAntennas;
            obj.PHY.MaximumDopplerShift             = MaximumDopplerShift;
            obj.Implementation.PowerDelayProfile    = PowerDelayProfile;
            obj.Implementation.timeCorrelation      = timeCorrelation;
            obj.Implementation.spatialCorrelation   = spatialCorrelation;
            obj.PHY.K                               = K;
            obj.PHY.delta                           = Delta;

            
            % calculate spatial correlation matrix according to TS36.101, Annex B
            if ~strcmp(obj.Implementation.spatialCorrelation,'none')
                if obj.Nr.txAntennas == 1
                    obj.Implementation.spatialCorrelationMatrixTX = 1;
                else
                    obj.Implementation.spatialCorrelationMatrixTX = toeplitz(correlationCoeffTx.^(linspace(0,1,obj.Nr.txAntennas).^2));
                end
                
                if obj.Nr.rxAntennas == 1
                    obj.Implementation.spatialCorrelationMatrixRX = 1;
                else
                    obj.Implementation.spatialCorrelationMatrixRX = toeplitz(correlationCoeffRx.^(linspace(0,1,obj.Nr.rxAntennas).^2));
                end
            end
            
            obj.PHY.dt = 1/obj.PHY.SamplingRate;
            
            if ischar(PowerDelayProfile)
                % For the new 3GPP 38.900 channel model, see http://www.3gpp.org/DynaReport/38900.htm
                if strcmp(PowerDelayProfile(1:3),'TDL')
                    Pos1 = strfind(PowerDelayProfile,'_');
                    Pos2 = strfind(PowerDelayProfile,'ns');
                    DesiredRMSdelaySpread = str2double(PowerDelayProfile(Pos1+1:Pos2-1))*10^-9;
                    PowerDelayProfile = PowerDelayProfile(1:5);
                end
                % See ftp://www.3gpp.org/tsg_ran/TSG_RAN/TSGR_16/Docs/PDF/RP-020376.pdf
                switch PowerDelayProfile
                    case 'Flat'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0;...                                           % Relative power in dB
                            0];                                             % Relative delay in ns
                    case 'AWGN'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0;...                                           % Relative power in dB
                            0];                                             % Relative delay in ns                        
                        if (obj.PHY.MaximumDopplerShift>0) && strcmp(PowerDelayProfile,'AWGN')
                            warning('The velocity is set to zero (AWGN channel!)');
                            obj.PHY.MaximumDopplerShift = 0;
                        end
                    case 'PedestrianA'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0,-9.7,-19.2,-22.8;...                          % Relative power in dB
                            0,110e-9,190e-9,410e-9];                        % Relative delay in ns
                    case 'PedestrianB'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0,-0.9,-4.9,-8,-7.8,-23.9;...                           
                            0,200e-9,800e-9,1200e-9,2300e-9,3700e-9];
                    case 'VehicularA'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0,-1,-9,-10,-15,-20;...                            
                            0,310e-9,710e-9,1090e-9,1730e-9,2510e-9];
                    case 'VehicularB'
                        PowerDelayProfileINdB_DelayIndex = [...
                            -2.5,0,-12.8,-10,-25.2,-16;...                            
                            0,300e-9,8900e-9,12900e-9,17100e-9,20000e-9];
                    case 'ExtendedPedestrianA'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0,-1,-2,-3,-8,-17.2,-20.8;
                            0,30e-9,70e-9,90e-9,110e-9,190e-9,410e-9];
                    case 'ExtendedVehicularA'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0,-1.5,-1.4,-3.6,-0.6,-9.1,-7,-12,-16.9; 
                            0,30e-9,150e-9,310e-9,370e-9,710e-9,1090e-9,1730e-9,2510e-9];
                    case {'ExtendedTypicalUrban', 'ETU'}
                        PowerDelayProfileINdB_DelayIndex = [...
                            -1,-1,-1,0,0,0,-3,-5,-7;
                            0,50e-9,120e-9,200e-9,230e-9,500e-9,1600e-9,2300e-9,5000e-9];
                    case {'HighSpeedTrain', 'HST'}
                        PowerDelayProfileINdB_DelayIndex = [...
                            -0.2;                       
                            0];
                    case 'TDL-A'
                        PowerDelayProfileINdB_DelayIndex(1,:) = ...
                            [-13.4,0,-2.2,-4,-6,-8.2,-9.9,-10.5,-7.5,-15.9,-6.6,-16.7,-12.4,-15.2,-10.8,-11.3,-12.7,-16.2,-18.3,-18.9,-16.6,-19.9,-29.7];
                        PowerDelayProfileINdB_DelayIndex(2,:) = ...
                            DesiredRMSdelaySpread*[0.0000,0.3819,0.4025,0.5868,0.4610,0.5375,0.6708,0.5750,0.7618,1.5375,1.8978,2.2242,2.1718,2.4942,2.5119,3.0582,4.0810,4.4579,4.5695,4.7966,5.0066,5.3043,9.6586];   
                    case 'TDL-B'
                        PowerDelayProfileINdB_DelayIndex(1,:) = ...
                            [0,-2.2,-4,-3.2,-9.8,-1.2,-3.4,-5.2,-7.6,-3,-8.9,-9,-4.8,-5.7,-7.5,-1.9,-7.6,-12.2,-9.8,-11.4,-14.9,-9.2,-11.3];
                        PowerDelayProfileINdB_DelayIndex(2,:) = ...
                            DesiredRMSdelaySpread*[0.0000,0.1072,0.2155,0.2095,0.2870,0.2986,0.3752,0.5055,0.3681,0.3697,0.5700,0.5283,1.1021,1.2756,1.5474,1.7842,2.0169,2.8294,3.0219,3.6187,4.1067,4.2790,4.7834];  
                    case 'TDL-C'
                        PowerDelayProfileINdB_DelayIndex(1,:) = ...
                            [-4.4,-1.2,-3.5,-5.2,-2.5,0,-2.2,-3.9,-7.4,-7.1,-10.7,-11.1,-5.1,-6.8,-8.7,-13.2,-13.9,-13.9,-15.8,-17.1,-16,-15.7,-21.6,-22.8];
                        PowerDelayProfileINdB_DelayIndex(2,:) = ...
                            DesiredRMSdelaySpread*[0,0.2099,0.2219,0.2329,0.2176,0.6366,0.6448,0.6560,0.6584,0.7935,0.8213,0.9336,1.2285,1.3083,2.1704,2.7105,4.2589,4.6003,5.4902,5.6077,6.3065,6.6374,7.0427,8.6523];  
                    case 'TDL-D'
                        PowerDelayProfileINdB_DelayIndex(1,:) = ...
                            [-0.2 -13.5 -18.8 -21 -22.8 -17.9 -20.1 -21.9 -22.9 -27.8 -23.6 -24.8 -30.0 -27.7];
                        PowerDelayProfileINdB_DelayIndex(2,:) = ...
                            DesiredRMSdelaySpread*[0 0 0.035 0.612 1.363 1.405 1.804 2.596 1.775 4.042 7.937 9.424 9.708 12.525];
                    case 'TDL-E'
                        PowerDelayProfileINdB_DelayIndex(1,:) = ...
                            [-0.03 -22.03 -15.8 -18.1 -19.8 -22.9 -22.4 -18.6 -20.8 -22.6 -22.3 -25.6 -20.2 -29.8 -29.2];
                        PowerDelayProfileINdB_DelayIndex(2,:) = ...
                            DesiredRMSdelaySpread*[0 0 0.5133 0.5440 0.5630 0.5440 0.7112 1.9092 1.9293 1.9589 2.6426 3.7136 5.4524 12.0034 20.6519];
                    case 'SUI-1'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0 -15 -20;
                            0 0.4e-6 0.9e-6];
                    case 'SUI-2'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0 -12 -15;
                            0 0.4e-6 1.1e-6];
                    case 'SUI-3'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0 -5 -10;
                            0 0.4e-6 0.9e-6];
                    case 'SUI-4'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0 -4 -8;
                            0 1.5e-6 4e-6];
                    case 'SUI-5'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0 -5 -10;
                            0 4e-6 10e-6];
                    case 'SUI-6'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0 -10 -14;
                            0 14e-6 20e-6];
                    case 'COST207-RA'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0 -4 -8 -12 -16 -20;
                            0 0.1e-6 0.2e-6 0.3e-6 0.4e-6 0.5e-6];
                    case 'COST207-TU'
                        PowerDelayProfileINdB_DelayIndex = [...
                            -3 0 -2 -6 -8 -10;
                            0 0.2e-6 0.5e-6 1.6e-6 2.3e-6 5.0e-6];
                    case 'COST207-BU'
                        PowerDelayProfileINdB_DelayIndex = [...
                            -7 -3 -1 0 -2 -6 -7 -1 -2 -7 -10 -15;
                            0 0.2e-6 0.4e-6 0.8e-6 1.6e-6 2.2e-6 3.2e-6 5.0e-6 6.0e-6 7.2e-6 8.2e-6 10.0e-6];
                    case 'COST207-HT'
                        PowerDelayProfileINdB_DelayIndex = [...
                            -10 -8 -6 -4 0 0 -4 -8 -9 -10 -12 -14;
                            0 0.2e-6 0.4e-6 0.6e-6 0.8e-6 2.0e-6 2.4e-6 15.0e-6 15.2e-6 15.8e-6 17.2e-6 20.0e-6];
                    otherwise
                        error('Power delay profile model not supported!');
                end
                IndexDelays = round(PowerDelayProfileINdB_DelayIndex(2,:)./obj.PHY.dt)+1;
                if  WarningIfSampleRateDoesNotMatch && (sum(abs(rem(PowerDelayProfileINdB_DelayIndex(2,:),obj.PHY.dt)))>0)
                    disp('Sampling rate does not match the predefined delays of the channel model!');
                    disp('Delay taps are changed according to (top = desired value; bottom = chosen value):');
                    disp([PowerDelayProfileINdB_DelayIndex(2,:);(IndexDelays-1)*obj.PHY.dt]);
                end
                PowerDelayProfileTemp = zeros(length(IndexDelays),max(IndexDelays));
                for i_delay = 1:length(IndexDelays)
                    PowerDelayProfileTemp(i_delay,IndexDelays(i_delay)) = 10.^(PowerDelayProfileINdB_DelayIndex(1,i_delay)/10);
                end
                obj.PHY.PowerDelayProfile = sum(PowerDelayProfileTemp,1);
                obj.PHY.DesiredPowerDelayProfiledB = PowerDelayProfileINdB_DelayIndex; 
            else
                obj.PHY.PowerDelayProfile = PowerDelayProfile;
                
                obj.PHY.DesiredPowerDelayProfiledB(1,:) = 10*log10(PowerDelayProfile);
                obj.PHY.DesiredPowerDelayProfiledB(2,:) = (0:size(obj.PHY.DesiredPowerDelayProfiledB,2)-1)*obj.PHY.dt;                
            end
            obj.Implementation.PowerDelayProfileNormalized = obj.PHY.PowerDelayProfile.'/sum(obj.PHY.PowerDelayProfile);
            
            obj.Implementation.IndexDelayTaps = find(obj.PHY.PowerDelayProfile');
            
            NrOnes = sum(obj.Nr.SamplesTotal:-1:(obj.Nr.SamplesTotal-length(obj.PHY.PowerDelayProfile)+1));
            obj.Implementation.MappingConvolutionMatrix = nan(NrOnes,2);
            for i_NrTaps = 1:length(obj.PHY.PowerDelayProfile)
                obj.Implementation.MappingConvolutionMatrix(sum(obj.Nr.SamplesTotal:-1:obj.Nr.SamplesTotal-i_NrTaps+2)+(1:obj.Nr.SamplesTotal-i_NrTaps+1),:) = [-i_NrTaps+1+(i_NrTaps:obj.Nr.SamplesTotal).' (i_NrTaps:obj.Nr.SamplesTotal).'];
            end
           
            % Faster mapping because zeros are not mapped
            obj.Implementation.MappingConvolutionMatrixFast = obj.Implementation.MappingConvolutionMatrix;
            for i_NrTaps = find(obj.PHY.PowerDelayProfile==0)
                obj.Implementation.MappingConvolutionMatrixFast(sum(obj.Nr.SamplesTotal:-1:obj.Nr.SamplesTotal-i_NrTaps+2)+(1:obj.Nr.SamplesTotal-i_NrTaps+1),:) = -1;
            end  
            obj.Implementation.MappingConvolutionMatrixFast(obj.Implementation.MappingConvolutionMatrixFast(:,1)==-1,:)=[];
            
            if obj.PHY.MaximumDopplerShift/(obj.PHY.SamplingRate/obj.Nr.SamplesTotal)<=0.5 && strcmp(DopplerModel(1:min(8,end)),'Discrete') && obj.PHY.MaximumDopplerShift>0
                warning('Discrete Doppler spectrum: The velocity is so low, that it is set to zero. Change to a "Jakes" or "Uniform" spectrum or increase the number of symbols in time!');
                obj.PHY.MaximumDopplerShift = 0;
            end            

            if obj.PHY.MaximumDopplerShift>0
                obj.Implementation.DopplerModel = DopplerModel;
                if strcmp(obj.Implementation.DopplerModel(1:min(8,end)),'Discrete') && ~obj.Implementation.timeCorrelation
                    obj.Implementation.UseDiscreteDopplerSpectrum = true;
                    % Calculate frequency spacing between discrete Doppler shifts
                    df = obj.PHY.SamplingRate/obj.Nr.SamplesTotal;

                    NrDopplerShifts = ceil(obj.PHY.MaximumDopplerShift/df);

                    IntervalPoints = df*((-NrDopplerShifts-1:NrDopplerShifts)+0.5);
                    IntervalPoints(IntervalPoints<=-obj.PHY.MaximumDopplerShift)=-obj.PHY.MaximumDopplerShift;
                    IntervalPoints(IntervalPoints>=obj.PHY.MaximumDopplerShift)=obj.PHY.MaximumDopplerShift;

                    switch obj.Implementation.DopplerModel
                        case 'Discrete-Jakes'
                            DiscreteDopplerSpectrum = asin(IntervalPoints(2:end)/obj.PHY.MaximumDopplerShift)-asin(IntervalPoints(1:end-1)/obj.PHY.MaximumDopplerShift); % calculate the integral of the continous Doppler spectrum over a certain interval
                        case 'Discrete-Uniform'
                            DiscreteDopplerSpectrum = IntervalPoints(2:end)-IntervalPoints(1:end-1); % calculate the integral of the continous Doppler spectrum over a certain interval
                        otherwise
                            error('Doppler spectrum not supported');
                    end    

                    % Normalize
                    DiscreteDopplerSpectrum = DiscreteDopplerSpectrum./sum(DiscreteDopplerSpectrum);

                    obj.Implementation.DiscreteDopplerSpectrum = repmat(DiscreteDopplerSpectrum.',1,length(obj.Implementation.IndexDelayTaps));            
                elseif strcmp(obj.Implementation.DopplerModel(1:min(8,end)),'Discrete') && obj.Implementation.timeCorrelation
                    error('Discrete Doppler spectrum: Correlated fading not supported with discrete Doppler spectrum!');
                else
                    obj.Nr.Paths = Paths;                    
                    obj.Implementation.UseDiscreteDopplerSpectrum = false;
                    
                    if obj.Implementation.timeCorrelation
                        % generate independent random variables for time correlated fading channel model
                        obj.Implementation.randomVariables.theta    = 2*pi*rand(length(obj.Implementation.IndexDelayTaps),1,obj.Nr.rxAntennas,obj.Nr.txAntennas)-pi;
                        obj.Implementation.randomVariables.phi      = 2*pi*rand(length(obj.Implementation.IndexDelayTaps),1,obj.Nr.Paths,obj.Nr.rxAntennas,obj.Nr.txAntennas)-pi;
                        obj.Implementation.randomVariables.psi      = 2*pi*rand(length(obj.Implementation.IndexDelayTaps),1,obj.Nr.Paths,obj.Nr.rxAntennas,obj.Nr.txAntennas)-pi;
                    end
                end
            end
  
            obj.NewRealization(0);
            
            obj.Implementation.CancelElementsConvolutionMatrix = ones(obj.Nr.SamplesTotal,length(obj.PHY.PowerDelayProfile))==1;
            for i_NrTaps = 1:length(obj.PHY.PowerDelayProfile)-1
                obj.Implementation.CancelElementsConvolutionMatrix(i_NrTaps,(end-length(obj.PHY.PowerDelayProfile)+1+i_NrTaps):end) = false;
            end
            obj.Implementation.CancelElementsConvolutionMatrixFast = obj.Implementation.CancelElementsConvolutionMatrix(:,obj.Implementation.IndexDelayTaps);
                 
        end
        
        function NewRealization( obj, frame )
            % obtain new realization of the channel

            if strcmp(obj.Implementation.PowerDelayProfile, 'AWGN')
                obj.ImpulseResponse(1,1,:,:) = ones( obj.Nr.rxAntennas, obj.Nr.txAntennas );
            else
                % time-variant channel
                if obj.PHY.MaximumDopplerShift>0
                    obj.ImpulseResponse = zeros( obj.Nr.SamplesTotal, size(obj.Implementation.PowerDelayProfileNormalized,1), obj.Nr.rxAntennas, obj.Nr.txAntennas  );
                    if obj.Implementation.UseDiscreteDopplerSpectrum % Faster calulation using IFFT. Underlying assumption: discrete Doppler shifts => valid only for a high number of samples/Doppler shift
                        for iTx = 1:obj.Nr.txAntennas
                            for iRx = 1:obj.Nr.rxAntennas
                                NrDopplerShifts = (size(obj.Implementation.DiscreteDopplerSpectrum,1)-1)/2;

                                GaussUncorr1 = obj.Nr.SamplesTotal/sqrt(2)*(randn(NrDopplerShifts+1,length(obj.Implementation.IndexDelayTaps))+1j*randn(NrDopplerShifts+1,length(obj.Implementation.IndexDelayTaps)));
                                GaussUncorr2 = obj.Nr.SamplesTotal/sqrt(2)*(randn(NrDopplerShifts,length(obj.Implementation.IndexDelayTaps))+1j*randn(NrDopplerShifts,length(obj.Implementation.IndexDelayTaps)));    
                     
                                GaussUncorr1 = GaussUncorr1.*repmat(sqrt(obj.Implementation.PowerDelayProfileNormalized(obj.Implementation.IndexDelayTaps).'),size(GaussUncorr1,1),1);
                                GaussUncorr2 = GaussUncorr2.*repmat(sqrt(obj.Implementation.PowerDelayProfileNormalized(obj.Implementation.IndexDelayTaps).'),size(GaussUncorr2,1),1);
                                
                                obj.ImpulseResponse(:,obj.Implementation.IndexDelayTaps,iRx,iTx) = ifft(...
                                     [sqrt(obj.Implementation.DiscreteDopplerSpectrum(NrDopplerShifts+1:end,:)).*GaussUncorr1;...
                                      zeros(obj.Nr.SamplesTotal-NrDopplerShifts*2-1,length(obj.Implementation.IndexDelayTaps));...
                                      sqrt(obj.Implementation.DiscreteDopplerSpectrum(1:NrDopplerShifts,:)).*GaussUncorr2]);    
                                                                
                  
                            end
                        end
                    else % This is the true Doppler spectum (continiuous and not discrete!)
                        if obj.Implementation.timeCorrelation
                            % correlated Rayleigh fading channel according to:
                            % Yahong Rosa Zheng et. al., "Simulation Models With Correct Statistical Properties for Rayleigh Fading Channels", IEEE Trans. on Comm., Vol. 51, No. 6, June, 2003
                            % with a slight modification as described in the appendix of:
                            % Thomas Zemen et. al., "Time-Variant Channel Estimation Using Discrete Prolate Spheroidal Sequences", IEEE Trans. on Signal Processing, Vol. 53, No. 9, Sep. 2005
                            switch obj.Implementation.DopplerModel
                                case 'Jakes'
                                    ImpulseResponseTemp = zeros(length(obj.PHY.PowerDelayProfile),obj.Nr.SamplesTotal,obj.Nr.rxAntennas,obj.Nr.txAntennas);

                                    t = (0:obj.PHY.dt:obj.PHY.dt*(obj.Nr.SamplesTotal-1))+obj.PHY.dt*obj.Nr.SamplesTotal*frame;
                                    alpha = (2*pi*repmat(reshape((1:obj.Nr.Paths),[1,1,obj.Nr.Paths]),[length(obj.Implementation.IndexDelayTaps),1,1,obj.Nr.rxAntennas,obj.Nr.txAntennas])-pi+permute(repmat(obj.Implementation.randomVariables.theta,[1,1,1,1,obj.Nr.Paths]),[1 2 5 3 4])) ./ (4*obj.Nr.Paths);
                                    Xc = cos(obj.Implementation.randomVariables.psi).*cos(2*pi*obj.PHY.MaximumDopplerShift .* t.*cos(alpha)+obj.Implementation.randomVariables.phi); % bsxfun implicit since matlab 2016b
                                    Xs = sin(obj.Implementation.randomVariables.psi).*cos(2*pi*obj.PHY.MaximumDopplerShift .* t.*cos(alpha)+obj.Implementation.randomVariables.phi);
                                    ImpulseResponseTemp(obj.Implementation.IndexDelayTaps,:,:,:) = 2/sqrt(2*obj.Nr.Paths) * sum(Xc+1i*Xs,3);

                                    obj.ImpulseResponse = permute(obj.Implementation.PowerDelayProfileNormalized .* ImpulseResponseTemp,[2 1 3 4]);
                                case 'Uniform'
                                    error('Uniform Doppler spectrum currently not implemented for correlated fading!');
                                otherwise
                                    error('Doppler spectrum not supported');
                            end
                        else
                            % uncorrelated fading
                            ImpulseResponseTemp = zeros(length(obj.PHY.PowerDelayProfile),obj.Nr.SamplesTotal);
                            for iTx = 1:obj.Nr.txAntennas
                                for iRx = 1:obj.Nr.rxAntennas
                                    switch obj.Implementation.DopplerModel
                                        case 'Jakes'
                                            DopplerShifts = cos(rand([length(obj.Implementation.IndexDelayTaps), 1, obj.Nr.Paths])*2*pi)*obj.PHY.MaximumDopplerShift;
                                        case 'Uniform'
                                            DopplerShifts = 2*(rand([length(obj.Implementation.IndexDelayTaps), 1, obj.Nr.Paths])-0.5)*obj.PHY.MaximumDopplerShift;
                                        otherwise
                                            error('Doppler spectrum not supported');
                                    end
                                    RandomPhase = rand([length(obj.Implementation.IndexDelayTaps) 1 obj.Nr.Paths]);
                                    t = 0:obj.PHY.dt:obj.PHY.dt*(obj.Nr.SamplesTotal-1);
                                    ImpulseResponseTemp(obj.Implementation.IndexDelayTaps,:) = sum(exp(1j*2*pi*(bsxfun(@plus,RandomPhase,bsxfun(@times,DopplerShifts,t)))),3)./sqrt(obj.Nr.Paths);

                                    obj.ImpulseResponse(:,:,iRx,iTx) = bsxfun(@times,sqrt(obj.Implementation.PowerDelayProfileNormalized), ImpulseResponseTemp).';
                                end
                            end
                        end
                    end
                else
                    % time-invariant channel
                    sigD2=1/(2*(1+obj.PHY.K));
                    V1=sqrt(obj.PHY.K*sigD2/2)*(sqrt(1+obj.PHY.delta)+sqrt(1-obj.PHY.delta));
                    V2=sqrt(obj.PHY.K*sigD2/2)*(sqrt(1+obj.PHY.delta)-sqrt(1-obj.PHY.delta));
                    obj.ImpulseResponse = zeros(1,size(obj.Implementation.PowerDelayProfileNormalized,1), obj.Nr.rxAntennas, obj.Nr.txAntennas);
                    for iTx = 1:obj.Nr.txAntennas
                        for iRx = 1:obj.Nr.rxAntennas
                            obj.ImpulseResponse(1,:,iRx,iTx) = ((sqrt(sigD2/length(obj.Implementation.PowerDelayProfileNormalized))).*(randn(size(obj.Implementation.PowerDelayProfileNormalized))+1j*randn(size(obj.Implementation.PowerDelayProfileNormalized)))).';
                            obj.ImpulseResponse(1,1,iRx,iTx) = obj.ImpulseResponse(1,1,iRx,iTx) + V1*exp(1j*2*pi*rand()) + V2*exp(1j*2*pi*rand());
                        end
                    end
                end
            end
            
            % introduce spatial correlation via Kronecker correlation model
            if ~strcmp(obj.Implementation.spatialCorrelation,'none')
                for iTime = 1:size(obj.ImpulseResponse,1)
                    for iDelay = 1:size(obj.ImpulseResponse,2)
                        obj.ImpulseResponse(iTime,iDelay,:,:) = sqrtm(obj.Implementation.spatialCorrelationMatrixRX) * reshape(obj.ImpulseResponse(iTime,iDelay,:,:),obj.Nr.rxAntennas, obj.Nr.txAntennas) * sqrtm(obj.Implementation.spatialCorrelationMatrixTX);
                    end
                end
            end
        end 
        
        
        function convolvedSignal = Convolution(obj, signal)
            % Convolution of the signal with the time variant impulse response
            
            N = size(signal,1);
            if obj.PHY.MaximumDopplerShift>0
                convolvedSignal = zeros(N, obj.Nr.rxAntennas);
                ConvolutionMatrix   = obj.GetConvolutionMatrix;
                for iTx = 1:obj.Nr.txAntennas
                    for iRx = 1:obj.Nr.rxAntennas
                        convolvedSignal(:,iRx) = convolvedSignal(:,iRx) + ConvolutionMatrix{iRx,iTx}(1:N,1:N)*signal(:,iTx);
                    end
                end
            else
                convolvedSignal = zeros(N + size(obj.Implementation.PowerDelayProfileNormalized,1)-1, obj.Nr.rxAntennas);
                for iTx = 1:obj.Nr.txAntennas
                    for iRx = 1:obj.Nr.rxAntennas
                        convolvedSignal(:,iRx) = convolvedSignal(:,iRx) + conv( signal(:,iTx), obj.ImpulseResponse(1,:,iRx,iTx) );
                    end
                end
                convolvedSignal = convolvedSignal(1:N,:);% remove tails
            end
        end
        
        function ConvolutionMatrix = GetConvolutionMatrix(obj)
            % returns the time-variant convolution matrix
            
            ConvolutionMatrix = cell( obj.Nr.rxAntennas, obj.Nr.txAntennas );
            if obj.PHY.MaximumDopplerShift>0
                for iTx = 1:obj.Nr.txAntennas
                    for iRx = 1:obj.Nr.rxAntennas
                        ImpulseResponseTemp = obj.ImpulseResponse(:,obj.Implementation.IndexDelayTaps,iRx,iTx);
                        ConvolutionMatrix{iRx,iTx} = sparse(obj.Implementation.MappingConvolutionMatrixFast(:,2),obj.Implementation.MappingConvolutionMatrixFast(:,1),ImpulseResponseTemp(obj.Implementation.CancelElementsConvolutionMatrixFast),obj.Nr.SamplesTotal,obj.Nr.SamplesTotal);
                    end
                end
            else
                for iTx = 1:obj.Nr.txAntennas
                    for iRx = 1:obj.Nr.rxAntennas
                        ImpulseResponseTemp = obj.ImpulseResponse(ones(obj.Nr.SamplesTotal,1),obj.Implementation.IndexDelayTaps,iRx,iTx);
                        ConvolutionMatrix{iRx,iTx} = sparse(obj.Implementation.MappingConvolutionMatrixFast(:,2),obj.Implementation.MappingConvolutionMatrixFast(:,1),ImpulseResponseTemp(obj.Implementation.CancelElementsConvolutionMatrixFast),obj.Nr.SamplesTotal,obj.Nr.SamplesTotal);
                    end
                end
            end
        end
        
        function ChannelTransferFunction = GetTransferFunction(obj,TimePos,FFTSize,ActiveSubcarrier)
            % Calculates the channel transfer function at certain time
            % positions. The first argument is a vector and gives the time
            % positions. The second argument represents the FFT size of the
            % underlying modulation format. The third argument (if
            % specified) determines which subcarriers are active (only 
            % these frequency positions are returned).


            ChannelTransferFunction = zeros( FFTSize, length(TimePos), obj.Nr.rxAntennas, obj.Nr.txAntennas );
            if obj.PHY.MaximumDopplerShift==0
                TimePos = ones(length(TimePos),1);
            end
            for iTx = 1:obj.Nr.txAntennas
                for iRx = 1:obj.Nr.rxAntennas
                    ImpulseTemp = [obj.ImpulseResponse(TimePos,:,iRx,iTx).';zeros(FFTSize-size(obj.ImpulseResponse,2),length(TimePos))];
                    ChannelTransferFunction(:,:,iRx,iTx) = fft(ImpulseTemp,[],1);
                end
            end  
            if exist('ActiveSubcarrier','var')
                ChannelTransferFunction = ChannelTransferFunction(ActiveSubcarrier,:,:,:);
            end
        end     
        
        function [TimeCorrelation,Time]=GetTimeCorrelation(obj)
            % returns the time-autocorrelation function of the channel.

            Time = -obj.PHY.dt*(obj.Nr.SamplesTotal-1):obj.PHY.dt:obj.PHY.dt*(obj.Nr.SamplesTotal-1);
            if obj.PHY.MaximumDopplerShift>0
                if obj.Implementation.UseDiscreteDopplerSpectrum
                    % Can be implemented similar as the frequency
                    % correlation. Needs to be done when time
                    warning('Exact solution only for a continuous Doppler spectrum.');
                end            
                Time = -obj.PHY.dt*(obj.Nr.SamplesTotal-1):obj.PHY.dt:obj.PHY.dt*(obj.Nr.SamplesTotal-1);
                switch obj.Implementation.DopplerModel
                    case 'Jakes'
                        TimeCorrelation = besselj(0,pi*2*obj.PHY.MaximumDopplerShift*Time);
                    case 'Uniform'
                        TimeCorrelation = sinc(2*obj.PHY.MaximumDopplerShift*Time);
                end
            else
                TimeCorrelation = ones(1,obj.Nr.SamplesTotal*2-1);
            end
        end
        
        function [FrequencyCorrelation,Frequency] = GetFrequencyCorrelation(obj)
            % returns the frequency-autocorrelation function of the channel.
            
            df = 1/(obj.Nr.SamplesTotal*obj.PHY.dt);
            
            FrequencyCorrelation = fft([obj.Implementation.PowerDelayProfileNormalized;zeros(obj.Nr.SamplesTotal-length(obj.Implementation.PowerDelayProfileNormalized),1)]);
            FrequencyCorrelation = circshift(FrequencyCorrelation,[ceil(obj.Nr.SamplesTotal/2) 1]);
            
            Frequency =((1:obj.Nr.SamplesTotal)-ceil(obj.Nr.SamplesTotal/2)-1)*df;
        end
        
        function MeanDelay = GetMeanDelay(obj)
            % returns the mean delay of the channel.
            
            Tau = (0:(length(obj.Implementation.PowerDelayProfileNormalized))-1)*obj.PHY.dt;
            MeanDelay = sum(Tau.*obj.Implementation.PowerDelayProfileNormalized.');
        end
        function RmsDelaySpread = GetRmsDelaySpread(obj)
            % returns the root mean squared doppler spread of the channel
            
            Tau = (0:(length(obj.Implementation.PowerDelayProfileNormalized))-1)*obj.PHY.dt;
            MeanDelay = obj.GetMeanDelay;
            RmsDelaySpread = sqrt(sum(Tau.^2.*obj.Implementation.PowerDelayProfileNormalized.')-MeanDelay.^2);
        end
        function CorrelationMatrix = GetCorrelationMatrix(obj)
            % returns the correlation matrix of the vectorized convolution 
            % matrix. Let y=H*x be the channel description with H beeing
            % the time-variant convolution matrix. Then, this method
            % returns the E{H(:)H(:)'}.

            IndexTimeCorrelation = bsxfun(@plus,obj.Nr.SamplesTotal+(0:(obj.Nr.SamplesTotal-1)).',(0:-1:-(obj.Nr.SamplesTotal-1)));
            TimeCorrelation =  obj.GetTimeCorrelation;
            CorrelationMatrixOrdered = (...
                repmat(TimeCorrelation(IndexTimeCorrelation),[size(obj.Implementation.PowerDelayProfileNormalized,1) 1]).*...
                kron(sparse(obj.Implementation.PowerDelayProfileNormalized),ones(size(IndexTimeCorrelation))));
            IndexCorrMatrixConv = bsxfun(@plus,(1:(obj.Nr.SamplesTotal+1):obj.Nr.SamplesTotal^2).',(0:1:(size(obj.ImpulseResponse,2)-1)));
            %           CorrelationMatrixOld = sparse(...
            %             repmat(reshape(IndexCorrMatrixConv,[],1),obj.Nr.SamplesTotal,1),...
            %             reshape(kron(IndexCorrMatrixConv.',ones(size(IndexCorrMatrixConv,1),1)),[],1),...
            %             CorrelationMatrixOrdered(:));
            % Memory efficient solution
            SplittingFactor = ceil(numel(IndexCorrMatrixConv)*obj.Nr.SamplesTotal*8/1e9);
            while mod(obj.Nr.SamplesTotal/SplittingFactor,1)
                SplittingFactor = SplittingFactor+1;
            end
            if SplittingFactor>1
                CorrelationMatrix = sparse(obj.Nr.SamplesTotal^2+size(IndexCorrMatrixConv,2)-1,obj.Nr.SamplesTotal^2+size(IndexCorrMatrixConv,2)-1);
                Row = repmat(reshape(IndexCorrMatrixConv,[],1),obj.Nr.SamplesTotal/SplittingFactor,1);
                for i_split = 1:SplittingFactor
                    Index1Temp = (1:obj.Nr.SamplesTotal/SplittingFactor)+(i_split-1)*obj.Nr.SamplesTotal/SplittingFactor;
                    Index2Temp = (1:numel(CorrelationMatrixOrdered)/SplittingFactor)+(i_split-1)*numel(CorrelationMatrixOrdered)/SplittingFactor;
                    
                    Column = reshape(kron(IndexCorrMatrixConv(Index1Temp,:).',ones(size(IndexCorrMatrixConv,1),1)),[],1);
                    CorrelationMatrix = CorrelationMatrix+...
                        sparse(Row,Column,CorrelationMatrixOrdered(Index2Temp),...
                        obj.Nr.SamplesTotal^2+size(IndexCorrMatrixConv,2)-1,...
                        obj.Nr.SamplesTotal^2+size(IndexCorrMatrixConv,2)-1);
                end
            else
                CorrelationMatrix = sparse(...
                    repmat(reshape(IndexCorrMatrixConv,[],1),obj.Nr.SamplesTotal,1),...
                    reshape(kron(IndexCorrMatrixConv.',ones(size(IndexCorrMatrixConv,1),1)),[],1),...
                    CorrelationMatrixOrdered(:));
            end
            CorrelationMatrix = CorrelationMatrix(1:obj.Nr.SamplesTotal^2,1:obj.Nr.SamplesTotal^2);
        end
        
        function [TimeCorrelation,Time] = PlotTimeCorrelation(obj,TimeSpacing)
            % Plot the time correlation. The input argument represents
            % the time-spacing and plots discrete points at multiples
            % of the time-spacing (only for presentational purpose)
            
            [TimeCorrelation,Time]=obj.GetTimeCorrelation;
            figure();
            plot(Time,abs(TimeCorrelation));
            ylabel('|Time Correlation|');
            xlabel('Time (s)');
            
            if not(exist('TimeSpacing','var'))
                TimePointSymbols = (-ceil(Time(end)/TimeSpacing):ceil(Time(end)/TimeSpacing))*TimeSpacing;
                hold on;
                Plot1 = stem(TimePointSymbols,abs(interp1(Time,TimeCorrelation,TimePointSymbols)),'black');
                legend(Plot1,{'TimeSpacing'});            
            end    
        end
        
        function [FrequencyCorrelation,Frequency] = PlotFrequencyCorrelation(obj,FrequencySpacing)
            % Plot the frequency correlation. The input argument represents
            % the time-spacing and plots discrete points at multiples
            % of the time-spacing (only for presentational purpose)
            
            [FrequencyCorrelation,Frequency]=obj.GetFrequencyCorrelation;
            figure();
            plot(Frequency,abs(FrequencyCorrelation));
            ylabel('|Frequency Correlation|');
            xlabel('Frequency (Hz)');
            
            if not(exist('FrequencySpacing','var'))
                FrequencyPointSymbols = (-ceil(Frequency(end)/FrequencySpacing):ceil(Frequency(end)/FrequencySpacing))*FrequencySpacing;
                hold on;
                Plot1 = stem(FrequencyPointSymbols,abs(interp1(Frequency,FrequencyCorrelation,FrequencyPointSymbols)),'black');
                legend(Plot1,{'FrequencySpacing'});
            end
        end
        function [PowerDelayProfile,Tau] = PlotPowerDelayProfile(obj)
            % Plot the desired power delay profile and the chosen power delay
            % profile, limited by the sampling rate
            
            Tau = (0:(length(obj.Implementation.PowerDelayProfileNormalized))-1)*obj.PHY.dt;
            PowerDelayProfile = obj.Implementation.PowerDelayProfileNormalized;
            if nargout==0
                DesiredTemp = 10.^(obj.PHY.DesiredPowerDelayProfiledB(1,:)/10);
                DesiredTemp = DesiredTemp/sum(DesiredTemp);
                RMSDelaySpread = obj.GetRmsDelaySpread;
                stem(obj.PHY.DesiredPowerDelayProfiledB(2,:)/1e-6,DesiredTemp,'-x red');
                hold on;
                stem(Tau/1e-6,PowerDelayProfile,'-o  blue');
                xlim([-obj.PHY.dt/1e-6/2 Tau(end)/1e-6+obj.PHY.dt/1e-6/2]);
                ylim([0 1]);
                xlabel('Delay [µs]');
                ylabel('Power Delay Profile');
                title(['RMS Delay Spread: ' num2str(round(RMSDelaySpread/1e-9)) 'ns']);
                legend({'Desired Delay Taps','Chosen Delay Taps (due to sampling)'});  
            end
        end
        
    end
end
