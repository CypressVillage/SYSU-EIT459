classdef SimulationParameters < handle
    % this class contains all simulation parameters, including user defined
    % parameters and constants
    
    properties
        simulation      % general simulation parameters
        topology        % network topology, describing which nodes are connected to each other
        schedule        % the time constant user allocation
        modulation      % modulation related parameters
        phy             % physical parameters
        coding          % channel coding parameters
        channel         % parameters related to the channel model
        constants       % physical constants
        feedback        % feedback related parameters
        layerMapping 
        
    end
    
    methods
        function obj = SimulationParameters( scenario )
            % simulation parameters class constructor
            
            % load parameters from a scenario file
            obj.loadScenario( scenario );
            
            % pre-defined parameters that should not be changed
            obj.setConstantParameters;         
            
            % check parameters
            obj.checkParameters();
            
            % calculate all dependent parameters
            obj.dependentParameters();
            
        end %function SimulationParameters
        
        function dependentParameters(obj)
            % this method calculates all parameters that are dependent on
            % other parameters
            
            % calculate thermal noise power (k*T*B)
            obj.phy.noisePower = obj.constants.BOLTZMANN * obj.phy.temperature * obj.modulation.samplingRate;
            
            % set number of symbols for data transmission to total number
            % of symbols minus number of symbols used as time guard (CP)
            obj.modulation.MCSymbols = obj.modulation.nSymbolsTotal - obj.modulation.nGuardSymbols;
            
            % calculate the total number of samples, assuming there is at least one base station
            tmpMod                              =    Modulation.Modulator(  'QAM', obj.modulation.mcsValues(obj.modulation.mcs(1)).modulationOrder, 'OFDM',...
                                                                            {'none'},...                                    % data symbol spreading type
                                                                            false,...                                       % noise power estimation
                                                                            1,...                                           % Number subcarriers
                                                                            obj.modulation.MCSymbols(1),...                 % Number OFDM Symbols
                                                                            obj.modulation.subcarrierSpacing(1),...         % Subcarrier spacing (Hz)
                                                                            obj.modulation.samplingRate,...                 % Sampling rate (Samples/s)
                                                                            0,...                                           % Intermediate frequency first subcarrier (Hz)
                                                                            false,...                                       % Transmit real valued signal
                                                                            obj.modulation.nGuardSymbols(1)/...             % Cyclic prefix length (s)
                                                                            (obj.modulation.subcarrierSpacing(1)*...
                                                                            obj.modulation.MCSymbols(1)), ...  
                                                                            0 ...                                           % Zero guard length (s)
                                                                            );
                          
%             obj.modulation.totalNumberOfSamples = tmpMod.WaveformObject.Nr.SamplesTotal;
            
            % get spatial correlation parameter according to TS36.101, Annex B
            switch obj.channel.spatialCorrelation
                case {'none','low'}
                    obj.channel.spatialCorrBS = 0;
                    obj.channel.spatialCorrUE = 0;
                case 'medium'
                    obj.channel.spatialCorrBS = 0.3;
                    obj.channel.spatialCorrUE = 0.9;
                case 'high'
                    obj.channel.spatialCorrBS = 0.9;
                    obj.channel.spatialCorrUE = 0.9;
                otherwise
                    error('Level of spatial correlation not supported!');
            end
            
           
        end
        
        function checkParameters(obj)
            % this method check whether the choice of parameters is
            % consistent and allowed
            
            % check nodes
            if ischar(obj.topology.nodes) 
                tempStr = strsplit(obj.topology.nodes, ',');
                
                if length(tempStr) == 1
                    error('Topology nodes input invalid!');
                else
                    nUE = 0;
                    nBS = 0;
                    for iNode = 1:length(tempStr)
                        switch tempStr{iNode}(1:2)
                            case {'UE'}
                                if str2double(tempStr{iNode}(3:end)) > 0
                                    nUE = nUE + 1;
                                else
                                    error('Topology nodes input invalid!');
                                end
                            case {'BS'}
                                if str2double(tempStr{iNode}(3:end)) > 0
                                    nBS = nBS + 1;
                                else
                                    error('Topology nodes input invalid!');
                                end
                            otherwise
                                error('Topology nodes input invalid!');
                        end
                    end
                end
                
            else
                error('Topology nodes input invalid format!');
            end
            obj.topology.nNodes = nUE + nBS;
            
            % Get number of primary links
            obj.topology.nLinks = length(strsplit(obj.topology.primaryLinks, ','));
            
            % Check if sweep parameter is correct
            % unpack struct from property
            variableName = strsplit(obj.simulation.sweepParam{1},'.');
            
            if isprop(obj, variableName{1})
                % extract struct if it exists
                tmpStruct = obj.(variableName{1});
            else
                error('sweep parameter is not a valid parameter');
            end %if isprop
            
            if isempty(tmpStruct.(variableName{2}))
                % error when variable does not exist
                error('sweep parameter is not a valid parameter');
            end
            
            % Identifiy whether the sweep parameter is per Node, per BS,
            % per UE or per Link.
            
            switch variableName{2}
                case {'txPowerBaseStation', 'nAntennasBaseStation'}
                    sweepParamType = 'per BS';
                case {'txPowerUser', 'nAntennasUser', 'userVelocity'}
                    sweepParamType = 'per UE';
                case {'pathloss','attenuation'}    
                    sweepParamType = 'per Link';
            end
            
            % </
            % Fix the length of the simulation parameters (FFW, makes sense
            % to merge it with the checks below). This is useful when we
            % want to apply the same parameter to all links or nodes and do
            % not want to repeat it for each one. Rather just enter it once
            % and the following code takes care of the repetition.
            
            % Per BS parameters
            if length(obj.simulation.txPowerBaseStation) == 1; obj.simulation.txPowerBaseStation = repmat(obj.simulation.txPowerBaseStation, 1, nBS); end
            if length(obj.simulation.nAntennasBaseStation) == 1; obj.simulation.nAntennasBaseStation = repmat(obj.simulation.nAntennasBaseStation, 1, nBS); end
            if length(obj.simulation.amplifierOBO) == 1; obj.simulation.amplifierOBO = repmat(obj.simulation.amplifierOBO, 1, nBS); end
            if length(obj.simulation.smoothnessFactor) == 1; obj.simulation.smoothnessFactor = repmat(obj.simulation.smoothnessFactor, 1, nBS); end
            
            % modulation
            if length(obj.modulation.waveform) == 1; obj.modulation.waveform = repmat(obj.modulation.waveform, 1, nBS); end
            if length(obj.modulation.numerOfSubcarriers) == 1; obj.modulation.numerOfSubcarriers = repmat(obj.modulation.numerOfSubcarriers, 1, nBS); end
            if length(obj.modulation.subcarrierSpacing) == 1; obj.modulation.subcarrierSpacing = repmat(obj.modulation.subcarrierSpacing, 1, nBS); end
            if length(obj.modulation.nSymbolsTotal) == 1; obj.modulation.nSymbolsTotal = repmat(obj.modulation.nSymbolsTotal, 1, nBS); end
            if length(obj.modulation.nGuardSymbols) == 1; obj.modulation.nGuardSymbols = repmat(obj.modulation.nGuardSymbols, 1, nBS); end

            % Per UE parameters
            if length(obj.simulation.txPowerUser) == 1; obj.simulation.txPowerUser = repmat(obj.simulation.txPowerUser, 1, nUE); end
            if length(obj.simulation.nAntennasUser) == 1; obj.simulation.nAntennasUser = repmat(obj.simulation.nAntennasUser, 1, nUE); end
            if length(obj.simulation.userVelocity) == 1; obj.simulation.userVelocity = repmat(obj.simulation.userVelocity, 1, nUE); end
            
            % Per Link parameters
            if length(obj.simulation.pathloss) == 1; obj.simulation.pathloss = repmat(obj.simulation.pathloss, 1, obj.topology.nLinks); end
            if length(obj.modulation.nStreams) == 1; obj.modulation.nStreams = repmat(obj.modulation.nStreams, 1, obj.topology.nLinks); end
            if length(obj.modulation.precodingMatrix) == 1
                precodingMatrix = obj.modulation.precodingMatrix;
                obj.modulation.precodingMatrix = cell(obj.topology.nLinks,1);
                obj.modulation.precodingMatrix(:) = precodingMatrix;            
            end
            if length(obj.modulation.mcs) == 1; obj.modulation.mcs = repmat(obj.modulation.mcs, 1, obj.topology.nLinks); end
            
            % coding
            if length(obj.coding.code) == 1; obj.coding.code = repmat(obj.coding.code , 1, obj.topology.nLinks); end
            if length(obj.coding.decoding) == 1; obj.coding.decoding = repmat(obj.coding.decoding , 1, obj.topology.nLinks); end
            if length(obj.coding.decodingIterations) == 1; obj.coding.decodingIterations = repmat(obj.coding.decodingIterations , 1, obj.topology.nLinks); end
            % />
            
            % check simulation parameters
            if length(obj.simulation.txPowerUser) ~= nUE
                error('There must be a tx power defined for each user!');
            end
            if length(obj.simulation.txPowerBaseStation) ~= nBS
                error('There must be a tx power defined for each base station!');
            end
            
            if length(obj.simulation.nAntennasUser) ~= nUE
                error('There must be a number of antennas defined for each user!');
            end
            
            if length(obj.simulation.nAntennasBaseStation) ~= nBS
                error('There must be a number of antennas defined for base station!');
            end
            
            if ~islogical(obj.simulation.downlinkNonlinearity)
                error('Invalid nonlinearity selection parameter! Input must be a logical true/false!')
            end
            
            % check modulation parameters
            if length(obj.modulation.waveform) ~= nBS
                error('There must be a waveform defined for each base station!');
            end
            
            if length(obj.modulation.subcarrierSpacing) ~= nBS
                error('There must be a subcarrier spacing defined for each base station!');
            end
            
            if length(obj.modulation.numerOfSubcarriers) ~= nBS
                error('There must be a number of subcarriers defined for each base station!');
            end
            
            if length(obj.modulation.nSymbolsTotal) ~= nBS
                error('There must be a total number of symbols defined for each base station!');
            end
            
            if length(obj.modulation.nGuardSymbols) ~= nBS
                error('There must be a number of guard symbols defined for each base station!');
            end
            
            if length(obj.simulation.userVelocity) ~= nUE
                error('There must be user velocity defined for each user!');
            end
            
            % check channel coding parameters
            if length(obj.coding.code) ~= obj.topology.nLinks
                error('There must be a code defined for each link!');
            end
            
            if length(obj.coding.decoding) ~= obj.topology.nLinks
                error('There must be a decoding algorithm  defined for each link!');
            end
            
            if length(obj.coding.decodingIterations) ~= obj.topology.nLinks
                error('There must be a decoding iterations defined for each link!');
            end
            
             % set the feedback according to the transmission mode
            switch obj.modulation.transmissionMode
                case 'CLSM'
                    if ~obj.feedback.enable 
                        warning('For the CLSM mode feedback is activated.');
                        obj.feedback.enable                    = true;
                    end
                    if ~obj.feedback.pmi
                        warning('For the CLSM mode the PMI, RI are activated');
                        obj.feedback.pmi                       = true;
                        obj.feedback.ri                        = false;
                    end
                    if ~obj.feedback.cqi
                        warning('For the CLSM mode the CQI is activated');
                        obj.feedback.cqi                       = true;
                    end
                    if obj.feedback.delay>0 && ~obj.channel.timeCorrelation
                        warning('For a temporal uncorrelated channel and a delay larger than 0 activating the feedback does not yield any improvements'); 
                    end
                                               
                case 'OLSM'
                    if ~obj.feedback.enable 
                        warning('For the OLSM mode feedback is activated.');
                        obj.feedback.enable                    = true;
                    end
                    if obj.feedback.pmi
                        warning('For the OLSM mode the PMI is deactivated');
                        obj.feedback.pmi                       = false;
                    end
                    if ~obj.feedback.ri
                        warning('For the OLSM mode the RI is activated');
                        obj.feedback.ri                       = true;
                    end
                    if ~obj.feedback.cqi
                        warning('For the OLSM mode the CQI is activated');
                        obj.feedback.cqi                       = true;
                    end
                    if obj.feedback.delay>0 && ~obj.channel.timeCorrelation
                        warning('For a temporal uncorrelated channel and a delay larger than 0 activating the feedback does not yield any improvements'); 
                    end

                case 'TxD'
                    if obj.feedback.enable || obj.feedback.pmi || obj.feedback.ri || obj.feedback.cqi
                        warning('For the TxD the feedback is deactivated');
                        obj.feedback.enable                    = false;
                        obj.feedback.pmi                       = false;
                        obj.feedback.ri                        = false;
                        obj.feedback.cqi                       = false;
                    end
                    % only one stream for TxD
                    if (sum(obj.modulation.nStreams>1)>0)
                        warning('For TxD the number of streams was set to 1');
                        obj.modulation.nStreams                = ones(size(obj.modulation.nStreams));
                    end 
                    if (sum(obj.simulation.nAntennasBaseStation==1)>0)
                        warning('TxD with one transmit antenna does not add any diversity');
                    end
                otherwise
                    if obj.feedback.enable && obj.feedback.delay>0 && ~obj.channel.timeCorrelation
                        warning('For a temporal uncorrelated channel and a delay larger than 0 activating the feedback does not yield any improvements'); 
                    end
                    % feedback setting are set as defined in the scenario file
            end
            
            % Check precoding matrix, number of streams and mcs
            if  (~obj.feedback.enable || ~obj.feedback.pmi) && ~ obj.feedback.ri &&  ~strcmp(obj.modulation.transmissionMode,'TxD') 
                if obj.topology.nLinks ~= length(obj.modulation.nStreams)
                    error('You have defined %d links but scStr.modulation.nStreams contains only %d values. A value must be assigned to each primary link.',obj.topology.nLinks ,size(obj.modulation.nStreams,2));
                end
                if obj.topology.nLinks ~= length(obj.modulation.precodingMatrix)
                    error('You have defined %d links but in scStr.modulation.precodingMatrix only %d precoding matrices were defined. A precoding matrix must be assigned to each primary link.',obj.topology.nLinks ,size(obj.modulation.precodingMatrix,2));
                end
            else
                precodingMatrix = cell(obj.topology.nLinks,1);
                precodingMatrix(:) = {obj.modulation.precodingMatrix{1}};
                obj.modulation.precodingMatrix = precodingMatrix;
                obj.modulation.nStreams = repmat(obj.modulation.nStreams(1), 1, obj.topology.nLinks);
            end
            
            if  (~obj.feedback.enable || ~obj.feedback.ri) &&  ~strcmp(obj.modulation.transmissionMode,'TxD') 
                if obj.topology.nLinks ~= length(obj.modulation.nStreams)
                    error('You have defined %d links but scStr.modulation.nStreams contains only %d values. A value must be assigned to each primary link.',obj.topology.nLinks ,size(obj.modulation.nStreams,2));
                end
                if obj.topology.nLinks ~= length(obj.modulation.precodingMatrix)
                    error('You have defined %d links but in scStr.modulation.precodingMatrix only %d precoding matrices were defined. A precoding matrix must be assigned to each primary link.',obj.topology.nLinks ,size(obj.modulation.precodingMatrix,2));
                end
            else
                precodingMatrix = cell(obj.topology.nLinks,1);
                precodingMatrix(:) = {obj.modulation.precodingMatrix{1}};
                obj.modulation.precodingMatrix = precodingMatrix;
                obj.modulation.nStreams = repmat(obj.modulation.nStreams(1), 1, obj.topology.nLinks);
            end          
            if ~obj.feedback.enable || ~obj.feedback.cqi 
                if obj.topology.nLinks ~= length(obj.modulation.mcs)
                    error('You have defined %d links but scStr.modulation.mcs contains only %d values. A value must be assigned to each primary link.',obj.topology.nLinks ,size(obj.modulation.mcs,2));
                end
            else
                if obj.topology.nLinks ~= length(obj.modulation.mcs)
                    obj.modulation.mcs = ones(1,obj.topology.nLinks);
                end
            end
            
            % check plotResultsFor
            if length(obj.simulation.plotResultsFor) ~=1 && length(obj.simulation.plotResultsFor) ~= obj.topology.nNodes
                error('The parameter ''plotResultsFor'' has to either have a single entry (i.e., applied globally), or has number of entries equal to the total number of nodes in the topology.');
            end
            for iP = 1:length(obj.simulation.plotResultsFor)
                if obj.simulation.plotResultsFor(iP) ~= 1 && obj.simulation.plotResultsFor(iP) ~= 0
                    error('The entries of the parameter ''plotResultsFor'' can only take values of 0 or 1.');
                end
            end
            
            % check applySweepingTo
            switch sweepParamType
                case 'per Node'
                    if length(obj.simulation.applySweepingTo) ~= 1 && length(obj.simulation.applySweepingTo) ~= obj.topology.nNodes
                        error('The parameter ''applySweepingTo'' has to either have a single entry (i.e., applied globally), or has number of entries equal to the total number of nodes in the topology.');
                    end
                case 'per BS'
                    if length(obj.simulation.applySweepingTo) ~= 1 && length(obj.simulation.applySweepingTo) ~= nBS
                        error('The parameter ''applySweepingTo'' has to either have a single entry (i.e., applied globally), or has number of entries equal to the total number of BSs in the topology.');
                    end
                case 'per UE'
                    if length(obj.simulation.applySweepingTo) ~= 1 && length(obj.simulation.applySweepingTo) ~= nUE
                        error('The parameter ''applySweepingTo'' has to either have a single entry (i.e., applied globally), or has number of entries equal to the total number of UEs in the topology.');
                    end
                case 'per Link'
                    if length(obj.simulation.applySweepingTo) ~= 1 && length(obj.simulation.applySweepingTo) ~= obj.topology.nLinks
                        error('The parameter ''applySweepingTo'' has to either have a single entry (i.e., applied globally), or has number of entries equal to the total number of Links in the topology.');
                    end
            end           
            for iP = 1:length(obj.simulation.applySweepingTo)
                if obj.simulation.applySweepingTo(iP) ~= 1 && obj.simulation.applySweepingTo(iP) ~= 0
                    error('The entries of the parameter ''applySweepingTo'' can only take values of 0 or 1.');
                end
            end
            
            % Downlink schedule check
            if obj.simulation.simulateDownlink
                if length(obj.schedule.fixedScheduleDL) ~= nBS
                    error('There must be a schedule specified for each base station!');
                else
                    for iBS = 1:nBS
                        if ~isempty(obj.schedule.fixedScheduleDL{iBS})
                            tempStr = strsplit(obj.schedule.fixedScheduleDL{iBS},',');
                            nSubcarriersTemp = 0;
                            for iNode = 1:length(tempStr)
                                tmpNodeStr = strsplit(tempStr{iNode},':');
                                if length(tmpNodeStr) == 2
                                    if strcmp(tmpNodeStr{1}(1:2), 'UE')
                                        if ~strcmp(tmpNodeStr{2}(1:2), 'UE')
                                            if str2double(tmpNodeStr{2}) > 0
                                                if mod(str2double(tmpNodeStr{2}),12)>0 && strcmp(obj.simulation.pilotPattern, 'LTE Downlink')
                                                    error('When the LTE Downlink pilot pattern is selected, the number of scheduled subcarriers must be a mutiple of 12!');
                                                end
                                                nSubcarriersTemp = nSubcarriersTemp + str2double(tmpNodeStr{2});
                                            else
                                                error('Schedule of base station %d is invalid!',iBS);
                                            end
                                        end
                                    elseif strcmp(tmpNodeStr{1}(1:4), 'none')
                                        if str2double(tmpNodeStr{2}) >= 0
                                            nSubcarriersTemp = nSubcarriersTemp + str2double(tmpNodeStr{2});
                                        else
                                            error('Schedule of base station %d is invalid!',iBS);
                                        end
                                    else
                                        error('Schedule of base station %d is invalid!',iBS);
                                    end
                                else
                                    error('Schedule of base station %d is invalid!',iBS);
                                end
                            end           
                            if nSubcarriersTemp ~= obj.modulation.numerOfSubcarriers(iBS)
                                error('Scheduled number of subcarriers in base station %d must add to the total number of %d subcarriers.', iBS, obj.modulation.numerOfSubcarriers(iBS));
                            end
                        end
                    end
                end
            end
             
            % Uplink schedule check
            if obj.simulation.simulateUplink
                if length(obj.schedule.fixedScheduleUL) ~= nBS
                    error('There must be a schedule specified for each base station!');
                else                
                    for iBS = 1:nBS
                        if ~isempty(obj.schedule.fixedScheduleUL{iBS})
                            tempStr = strsplit(obj.schedule.fixedScheduleUL{iBS},',');
                            nSubcarriersTemp = 0;
                            for iNode = 1:length(tempStr)
                                tmpNodeStr = strsplit(tempStr{iNode},':');
                                if length(tmpNodeStr) == 2
                                    if strcmp(tmpNodeStr{1}(1:2), 'UE')
                                        if ~strcmp(tmpNodeStr{2}(1:2), 'UE')
                                            if str2double(tmpNodeStr{2}) > 0
                                                if mod(str2double(tmpNodeStr{2}),12)>0 && strcmp(obj.simulation.pilotPattern, 'LTE Downlink')
                                                    error('When the LTE Downlink pilot pattern is selected, the number of scheduled subcarriers must be a mutiple of 12!');
                                                end
                                                nSubcarriersTemp = nSubcarriersTemp + str2double(tmpNodeStr{2});
                                            else
                                                error('Schedule of base station %d is invalid!',iBS);
                                            end
                                        end
                                    elseif strcmp(tmpNodeStr{1}(1:4), 'none')
                                        if str2double(tmpNodeStr{2}) >= 0
                                            nSubcarriersTemp = nSubcarriersTemp + str2double(tmpNodeStr{2});
                                        else
                                            error('Schedule of base station %d is invalid!',iBS);
                                        end
                                    else
                                        error('Schedule of base station %d is invalid!',iBS);
                                    end
                                else
                                    error('Schedule of base station %d is invalid!',iBS);
                                end
                            end           
                            if nSubcarriersTemp ~= obj.modulation.numerOfSubcarriers(iBS)
                                error('Scheduled number of subcarriers in base station %d must add to the total number of %d subcarriers.', iBS, obj.modulation.numerOfSubcarriers(iBS));
                            end
                        end
                    end
                end
            end            
            
            % check relation of modulation parameters between cells
            % for FBMC the number of real valued symbols is twice as high compared to a non FBMC cell with complex valued symbols
            fbmc_index = zeros(nBS,1);
            for iBS = 1:nBS
                fbmc_index(iBS,1) = double(strcmp(obj.modulation.waveform(iBS),'FBMC'))+1;
            end
            
            for iBS = 2:nBS
                tmpFactor = obj.modulation.subcarrierSpacing(1) / obj.modulation.subcarrierSpacing(iBS);
                       
                if( obj.modulation.nSymbolsTotal(iBS) / fbmc_index(iBS) * tmpFactor ~= obj.modulation.nSymbolsTotal(1) / fbmc_index(1) )
                    error('Invalid ratio of subcarrier spacing, number of subcarriers and number of total symbols across cells.');
                end
            end
            
            % check and set sampling rate           
            if ischar(obj.modulation.samplingRate)
                if strcmp(obj.modulation.samplingRate,'Automatic')
                    % automatically calculate sampling rate     
                    % choose sampling frequency such that the guard length is an integer multiple of the sampling time
                    
                    % get the greatest common divisor of guard times of all cells
                    [~, guardLengthGCD_numerator, guardLengthGCD_denominator] = Utils.vectorGCD_rational(obj.modulation.nGuardSymbols, (obj.modulation.nSymbolsTotal - obj.modulation.nGuardSymbols) .* obj.modulation.subcarrierSpacing);
                    
                    % make sure that the minimum oversampling factor is according to the parameter settings
                    minimumSamplingRate         = max(obj.modulation.subcarrierSpacing) * max(obj.modulation.numerOfSubcarriers);
                    appliedOversamplingFactor   = max(2, guardLengthGCD_denominator / gcd(guardLengthGCD_numerator * minimumSamplingRate, guardLengthGCD_denominator)); % min oversampling should be 2
                    obj.modulation.samplingRate = appliedOversamplingFactor * minimumSamplingRate;

                    fprintf('The sampling rate was set to %f MHz.\n', obj.modulation.samplingRate/1e6);
                else
                    error('Sampling rate statement unknown!');
                end
            elseif isnumeric(obj.modulation.samplingRate)
                % check is there is only one sampling rate
                if length(obj.modulation.samplingRate) > 1
                    error('There must be a single sampling rate for all cells!');
                end
                % check if sampling rate is multiple of 1/guard time
                if any(mod(obj.modulation.samplingRate .* obj.modulation.nGuardSymbols ./ ((obj.modulation.nSymbolsTotal - obj.modulation.nGuardSymbols) .* obj.modulation.subcarrierSpacing), 1)>0)
                    warning('Guard time is not integer multiple of the sampling time!');
                end
            else
                error('Sampling rate input format invalid!');
            end
            
            % check channel
            
            if strcmp(obj.channel.powerDelayProfile, 'AWGN') && (any(obj.simulation.nAntennasUser > 1) || any(obj.simulation.nAntennasBaseStation > 1))
                error('AWGN channel not defined for MIMO transmission!');
            end
            
            if obj.channel.timeCorrelation && ~strcmp(obj.channel.dopplerModel,'Jakes')
                error('Time correlated fading only supported with the ''Jakes'' Doppler spectrum.');
            end
            
            if ~ischar(obj.channel.spatialCorrelation)
                error('Spatial channel correlation parameter must be string!');
            end
            
            if ~strcmp(obj.channel.spatialCorrelation,'none') && all(obj.simulation.nAntennasUser==1) && all(obj.simulation.nAntennasBaseStation==1)
                warning('Spatial channel correlation is set to %s while all nodes have only a single antenna. Spatial correlation has no effect in a SISO simulation!',obj.channel.spatialCorrelation);
            end
            
            if ~all(~obj.simulation.userVelocity) && obj.channel.K~=0
                warning('In the mixcure case of moving and stationary users, Rayleigh fading will be applied to the stationary ones!')
                obj.channel.K = 0;
            end
            
            if (obj.channel.delta < 0) || (obj.channel.delta >1)
                error('TWDP parameter delta must be between 0 and 1!');
            end
            
            if (obj.channel.K < 0)
                error('TWDP parameter K must be positive!');
            end
            
            % check if sampling rate fits channel model
            % the check is implemented in the channel class -> just
            % generate a dummy object
            Channel.FastFading( obj.modulation.samplingRate, ...        % sampling rate
                                obj.channel.powerDelayProfile, ...      % channel model
                                1, ...                                  % total number of samples
                                0, ...                                  % doppler frequency
                                obj.channel.dopplerModel,...            % doppler model
                                obj.channel.nPaths,...                  % number of paths
                                false,...                               % time correlated fading
                                'none',...                              % spatial correlation
                                0,...                                   % spatial correlation coeff TX
                                0,...                                   % spatial correlation coeff RX
                                1,...                                   % number of transmit antennas
                                1,...                                   % number of receive antennas
                                true,...                                % show checks
                                0,...                                   % K = 0
                                1 ...                                   % delta = 1
                                );


             % check feedback object
             if obj.feedback.enable
                 if (isfield(obj.feedback,'delay') && isfield(obj.feedback.averager,'Type'))
                     if obj.feedback.delay<0
                         error('Feedback delay can not be negative!');
                     end
                     if (~strcmp(obj.feedback.averager.Type,'eesm') && ~strcmp(obj.feedback.averager.Type,'miesm'))
                          error('Only EESM and MIESM averager are supported!');
                     end
                 else
                     error('To enable to feedback a delay and an averager type are necessary!');
                 end
             end
            
             if ~isfield(obj.simulation,'saveData')
                obj.simulation.saveData = false;
             end
        end
        
        function loadScenario(obj, scenario)
            % this method loads parameters from a scenario file
            
            % check input
            if ~ischar(scenario)
                error('Scenario must be string!');
            end
            if exist(['Scenarios/', scenario, '.m'], 'file') == 2
                run(['Scenarios/', scenario]);
            else
                error('Scenario does not exist!');
            end
            
            % copy scenario parameters to object
            obj.simulation      = scStr.simulation;
            obj.topology        = scStr.topology;
            obj.schedule        = scStr.schedule;
            obj.channel         = scStr.channel;
            obj.modulation      = scStr.modulation;
            obj.coding          = scStr.coding;
            obj.layerMapping    = scStr.layerMapping;
            if isfield(scStr,'feedback')
                obj.feedback=scStr.feedback;
            else 
                obj.feedback.enable=0;
            end
            
            
        end
        
        function UpdateSweepValue(obj,iSweep)
            % this method overwirtes a specific object property with the
            % current sweep variable value
            
            variableName = strsplit(obj.simulation.sweepParam{1},'.');
            if obj.simulation.applySweepingTo == 1
                obj.simulation.applySweepingTo = ones(size(obj.(variableName{1}).(variableName{2})));
            elseif sum(obj.simulation.applySweepingTo) == 0 || isempty(obj.simulation.applySweepingTo)
                error('Atleast one element in obj.simulation.applySweepingTo has to be set to 1.');
            end
            obj.(variableName{1}).(variableName{2})(logical(obj.simulation.applySweepingTo)) = obj.simulation.sweepValue(iSweep);
        end
        
        function setConstantParameters( obj )
            % defines constant and pre-defined parameters that should not
            % be changed
            
            % define physical constants in SI units
            obj.constants.SPEED_OF_LIGHT            = 299792458;
            obj.constants.BOLTZMANN                 = 1.38064852e-23;
            
            % set parameters that are assumed to be constant
            obj.phy.temperature                     = 293;
            switch obj.modulation.cqiTable 
                case 0
                    % TS 38.214 Table 5.2.2.1-2 (QPSK to 64 QAM)
                    % pre-defined MCS
                    obj.modulation.mcsValues(1).cqi = 1;
                    obj.modulation.mcsValues(1).alphabet = 'QPSK';
                    obj.modulation.mcsValues(1).modulationOrder = 4;
                    obj.modulation.mcsValues(1).codingRateTimes1024 = 78;
                    obj.modulation.mcsValues(1).efficiency = 0.1523;

                    obj.modulation.mcsValues(2).cqi = 2;
                    obj.modulation.mcsValues(2).alphabet = 'QPSK';
                    obj.modulation.mcsValues(2).modulationOrder = 4;
                    obj.modulation.mcsValues(2).codingRateTimes1024 = 120;
                    obj.modulation.mcsValues(2).efficiency = 0.2344;

                    obj.modulation.mcsValues(3).cqi = 3;
                    obj.modulation.mcsValues(3).alphabet = 'QPSK';
                    obj.modulation.mcsValues(3).modulationOrder = 4;
                    obj.modulation.mcsValues(3).codingRateTimes1024 = 193;
                    obj.modulation.mcsValues(3).efficiency = 0.3770;

                    obj.modulation.mcsValues(4).cqi = 4;
                    obj.modulation.mcsValues(4).alphabet = 'QPSK';
                    obj.modulation.mcsValues(4).modulationOrder = 4;
                    obj.modulation.mcsValues(4).codingRateTimes1024 = 308;
                    obj.modulation.mcsValues(4).efficiency = 0.6016;

                    obj.modulation.mcsValues(5).cqi = 5;
                    obj.modulation.mcsValues(5).alphabet = 'QPSK';
                    obj.modulation.mcsValues(5).modulationOrder = 4;
                    obj.modulation.mcsValues(5).codingRateTimes1024 = 449;
                    obj.modulation.mcsValues(5).efficiency = 0.8770;

                    obj.modulation.mcsValues(6).cqi = 6;
                    obj.modulation.mcsValues(6).alphabet = 'QPSK';
                    obj.modulation.mcsValues(6).modulationOrder = 4;
                    obj.modulation.mcsValues(6).codingRateTimes1024 = 602;
                    obj.modulation.mcsValues(6).efficiency = 1.1758;

                    obj.modulation.mcsValues(7).cqi = 7;
                    obj.modulation.mcsValues(7).alphabet = '16QAM';
                    obj.modulation.mcsValues(7).modulationOrder = 16;
                    obj.modulation.mcsValues(7).codingRateTimes1024 = 378;
                    obj.modulation.mcsValues(7).efficiency = 1.4766;

                    obj.modulation.mcsValues(8).cqi = 8;
                    obj.modulation.mcsValues(8).alphabet = '16QAM';
                    obj.modulation.mcsValues(8).modulationOrder = 16;
                    obj.modulation.mcsValues(8).codingRateTimes1024 = 490;
                    obj.modulation.mcsValues(8).efficiency = 1.9141;

                    obj.modulation.mcsValues(9).cqi = 9;
                    obj.modulation.mcsValues(9).alphabet = '16QAM';
                    obj.modulation.mcsValues(9).modulationOrder = 16;
                    obj.modulation.mcsValues(9).codingRateTimes1024 = 616;
                    obj.modulation.mcsValues(9).efficiency = 2.4063;

                    obj.modulation.mcsValues(10).cqi = 10;
                    obj.modulation.mcsValues(10).alphabet = '64QAM';
                    obj.modulation.mcsValues(10).modulationOrder = 64;
                    obj.modulation.mcsValues(10).codingRateTimes1024 = 466;
                    obj.modulation.mcsValues(10).efficiency = 2.7305;

                    obj.modulation.mcsValues(11).cqi = 11;
                    obj.modulation.mcsValues(11).alphabet = '64QAM';
                    obj.modulation.mcsValues(11).modulationOrder = 64;
                    obj.modulation.mcsValues(11).codingRateTimes1024 = 567;
                    obj.modulation.mcsValues(11).efficiency = 3.3223;

                    obj.modulation.mcsValues(12).cqi = 12;
                    obj.modulation.mcsValues(12).alphabet = '64QAM';
                    obj.modulation.mcsValues(12).modulationOrder =64;
                    obj.modulation.mcsValues(12).codingRateTimes1024 = 666;
                    obj.modulation.mcsValues(12).efficiency = 3.9023;

                    obj.modulation.mcsValues(13).cqi = 13;
                    obj.modulation.mcsValues(13).alphabet = '64QAM';
                    obj.modulation.mcsValues(13).modulationOrder = 64;
                    obj.modulation.mcsValues(13).codingRateTimes1024 = 772;
                    obj.modulation.mcsValues(13).efficiency = 4.5234;

                    obj.modulation.mcsValues(14).cqi = 14;
                    obj.modulation.mcsValues(14).alphabet = '64QAM';
                    obj.modulation.mcsValues(14).modulationOrder = 64;
                    obj.modulation.mcsValues(14).codingRateTimes1024 = 873;
                    obj.modulation.mcsValues(14).efficiency = 5.1152;

                    obj.modulation.mcsValues(15).cqi = 15;
                    obj.modulation.mcsValues(15).alphabet = '64QAM';
                    obj.modulation.mcsValues(15).modulationOrder = 64;
                    obj.modulation.mcsValues(15).codingRateTimes1024 = 948;
                    obj.modulation.mcsValues(15).efficiency = 5.5547;

                    obj.feedback.averager.EESMbetas     = [5.01,5.01,0.84,1.67,1.61,1.64,3.87,5.06,6.4,12.59,17.59,23.33,29.45,33.05,35.41];
                    obj.feedback.averager.MIESMbetas    = [3.07,4.41,0.6,1.16,1.06,1.06,0.87,1.01,1.04,1.03,1.11,1.01,1.07,1,1.05];
                    obj.feedback.CqiMappingTable        = [-500;-5.6;-3.5;-1.5;0.7;2.6;4.6;6.5;8.5;10.5;12.5;14.3;16.1;18.1;20.3;22]; 

               
                case 1
                    % 1: TS 38.214 Table 5.2.2.1-3 (QPSK to 256 QAM)
                    obj.modulation.mcsValues(1).cqi = 1;
                    obj.modulation.mcsValues(1).alphabet = 'QPSK';
                    obj.modulation.mcsValues(1).modulationOrder = 4;
                    obj.modulation.mcsValues(1).codingRateTimes1024 = 78;
                    obj.modulation.mcsValues(1).efficiency = 0.1523;

                    obj.modulation.mcsValues(2).cqi = 2;
                    obj.modulation.mcsValues(2).alphabet = 'QPSK';
                    obj.modulation.mcsValues(2).modulationOrder = 4;
                    obj.modulation.mcsValues(2).codingRateTimes1024 = 193;
                    obj.modulation.mcsValues(2).efficiency = 0.3770;

                    obj.modulation.mcsValues(3).cqi = 3;
                    obj.modulation.mcsValues(3).alphabet = 'QPSK';
                    obj.modulation.mcsValues(3).modulationOrder = 4;
                    obj.modulation.mcsValues(3).codingRateTimes1024 = 449;
                    obj.modulation.mcsValues(3).efficiency = 0.8770;

                    obj.modulation.mcsValues(4).cqi = 4;
                    obj.modulation.mcsValues(4).alphabet = '16QAM';
                    obj.modulation.mcsValues(4).modulationOrder = 16;
                    obj.modulation.mcsValues(4).codingRateTimes1024 = 378;
                    obj.modulation.mcsValues(4).efficiency = 1.4766;

                    obj.modulation.mcsValues(5).cqi = 5;
                    obj.modulation.mcsValues(5).alphabet = '16QAM';
                    obj.modulation.mcsValues(5).modulationOrder = 16;
                    obj.modulation.mcsValues(5).codingRateTimes1024 = 490;
                    obj.modulation.mcsValues(5).efficiency = 1.9141;

                    obj.modulation.mcsValues(6).cqi = 6;
                    obj.modulation.mcsValues(6).alphabet = '16QAM';
                    obj.modulation.mcsValues(6).modulationOrder = 16;
                    obj.modulation.mcsValues(6).codingRateTimes1024 = 616;
                    obj.modulation.mcsValues(6).efficiency = 2.4063;

                    obj.modulation.mcsValues(7).cqi = 7;
                    obj.modulation.mcsValues(7).alphabet = '64QAM';
                    obj.modulation.mcsValues(7).modulationOrder = 64;
                    obj.modulation.mcsValues(7).codingRateTimes1024 = 466;
                    obj.modulation.mcsValues(7).efficiency = 2.7305;

                    obj.modulation.mcsValues(8).cqi = 8;
                    obj.modulation.mcsValues(8).alphabet = '64QAM';
                    obj.modulation.mcsValues(8).modulationOrder = 64;
                    obj.modulation.mcsValues(8).codingRateTimes1024 = 567;
                    obj.modulation.mcsValues(8).efficiency = 3.3223;

                    obj.modulation.mcsValues(9).cqi = 9;
                    obj.modulation.mcsValues(9).alphabet = '64QAM';
                    obj.modulation.mcsValues(9).modulationOrder = 64;
                    obj.modulation.mcsValues(9).codingRateTimes1024 = 666;
                    obj.modulation.mcsValues(9).efficiency = 3.9023;

                    obj.modulation.mcsValues(10).cqi = 10;
                    obj.modulation.mcsValues(10).alphabet = '64QAM';
                    obj.modulation.mcsValues(10).modulationOrder = 64;
                    obj.modulation.mcsValues(10).codingRateTimes1024 = 772;
                    obj.modulation.mcsValues(10).efficiency = 4.5234;

                    obj.modulation.mcsValues(11).cqi = 11;
                    obj.modulation.mcsValues(11).alphabet = '64QAM';
                    obj.modulation.mcsValues(11).modulationOrder = 64;
                    obj.modulation.mcsValues(11).codingRateTimes1024 = 873;
                    obj.modulation.mcsValues(11).efficiency = 5.1152;

                    obj.modulation.mcsValues(12).cqi = 12;
                    obj.modulation.mcsValues(12).alphabet = '256QAM';
                    obj.modulation.mcsValues(12).modulationOrder = 256;
                    obj.modulation.mcsValues(12).codingRateTimes1024 = 711;
                    obj.modulation.mcsValues(12).efficiency = 5.5547;

                    obj.modulation.mcsValues(13).cqi = 13;
                    obj.modulation.mcsValues(13).alphabet = '256QAM';
                    obj.modulation.mcsValues(13).modulationOrder = 256;
                    obj.modulation.mcsValues(13).codingRateTimes1024 = 797;
                    obj.modulation.mcsValues(13).efficiency = 6.2266;

                    obj.modulation.mcsValues(14).cqi = 14;
                    obj.modulation.mcsValues(14).alphabet = '256QAM';
                    obj.modulation.mcsValues(14).modulationOrder = 256;
                    obj.modulation.mcsValues(14).codingRateTimes1024 = 885;
                    obj.modulation.mcsValues(14).efficiency = 6.9141;

                    obj.modulation.mcsValues(15).cqi = 15;
                    obj.modulation.mcsValues(15).alphabet = '256QAM';
                    obj.modulation.mcsValues(15).modulationOrder = 256;
                    obj.modulation.mcsValues(15).codingRateTimes1024 = 948;
                    obj.modulation.mcsValues(15).efficiency = 7.4063;

                    obj.feedback.averager.EESMbetas     = [5.01,0.84,1.61,3.87,5.06,6.4,12.59,17.59,23.33,29.45,33.05,100.3,105.89,130,143.8];
                    obj.feedback.averager.MIESMbetas    = [3.07,0.6,1.06,0.87,1.01,1.04,1.03,1.11,1.01,1.07,1,1.19,0.92,0.97,1.12];
                    obj.feedback.CqiMappingTable        = [-500;-3.5;0.7;4.6;6.5;8.5;10.5;12.5;14.3;16.1;18.1;20.3;21.2;23.8;25.6;27];

                otherwise
                    error('CQI table not defined');
            end
            
            obj.feedback.averager.MCS_values       = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];         % Modulation a Coding Schemes used for the transmission

            % Layer Mapping
            if strcmp(obj.layerMapping.mode,'LTE')
                obj.layerMapping.table.Uplink = {1;[1 1];[1 2];[2 2]}; 
                obj.layerMapping.table.Downlink = {1;[1 1];[1 2];[2 2];[2 3];[3 3];[3 4];[4 4]};
            elseif strcmp(obj.layerMapping.mode,'5G')
                obj.layerMapping.table.Uplink = {1;2;3;4}; 
                obj.layerMapping.table.Downlink = {1;2;3;4;[2 3];[3 3];[3 4];[4 4]}; 
            else
                error('Layer mapping mode unkown!');
            end
        end %setConstantParameters
        
        function Links = initializeLinks(obj, Links, BS, UE)
            % this function  initializes all links according to some
            % scenario setting
            L = length(Links);
            
            nBS = length( BS );
            
            % get schedule for each BS downlink
            if obj.simulation.simulateDownlink
                linkScheduleDL = cell(nBS, 1);
                for iBS = 1:nBS
                    if ~isempty(obj.schedule.fixedScheduleDL{iBS})
                        split = strsplit( obj.schedule.fixedScheduleDL{iBS}, ',' );
                        nFrequBlocks = length( split );
                        linkScheduleDL{iBS} = [];
                        for iBlock = 1:nFrequBlocks
                            split2 = strsplit( split{iBlock}, ':' );
                            if strcmp(split2{1}(1:2),'UE')
                                if ~strcmp(split2{2}(1:2),'UE')
                                    tmpSchedule = UE{str2double(split2{1}(3:end))}.ID * ones(str2double(split2{2}),1);
                                else
                                    Links{BS{iBS}.ID, UE{str2double(split2{2}(3:end))}.ID}.MUSTFarUE = UE{str2double(split2{1}(3:end))}.ID;
                                    Links{BS{iBS}.ID, UE{str2double(split2{1}(3:end))}.ID}.MUSTNearUE = UE{str2double(split2{2}(3:end))}.ID;
                                    Links{BS{iBS}.ID, UE{str2double(split2{2}(3:end))}.ID}.IsMUSTNearUE = 1;
                                    Links{BS{iBS}.ID, UE{str2double(split2{1}(3:end))}.ID}.IsMUSTFarUE = 1;          
                                    tmpSchedule = [];
                                end
                            elseif strcmp(split2{1}(1:4),'none')
                                tmpSchedule = zeros(str2double(split2{2}),1);
                            else
                                error('Schedule input invalid!');
                            end
                            linkScheduleDL{iBS} = [ linkScheduleDL{iBS}; tmpSchedule ];
                        end
                        linkScheduleDL{iBS} = repmat( linkScheduleDL{iBS}, 1, obj.modulation.MCSymbols(iBS) );
                    end
                end
            end
            
            % get schedule for each BS Uplink
            if obj.simulation.simulateUplink
                linkScheduleUL = cell(nBS, 1);
                for iBS = 1:nBS
                    if ~isempty(obj.schedule.fixedScheduleUL{iBS})
                        split = strsplit( obj.schedule.fixedScheduleUL{iBS}, ',' );
                        nFrequBlocks = length( split );
                        linkScheduleUL{iBS} = [];
                        for iBlock = 1:nFrequBlocks
                            split2 = strsplit( split{iBlock}, ':' );
                            if strcmp(split2{1}(1:2),'UE')
                                tmpSchedule = UE{str2double(split2{1}(3:end))}.ID * ones(str2double(split2{2}),1);
                            elseif strcmp(split2{1}(1:4),'none')
                                tmpSchedule = zeros(str2double(split2{2}),1);
                            else
                                error('Schedule input invalid!');
                            end
                            linkScheduleUL{iBS} = [ linkScheduleUL{iBS}; tmpSchedule ];
                        end
                        linkScheduleUL{iBS} = repmat( linkScheduleUL{iBS}, 1, obj.modulation.MCSymbols(iBS) );
                    end
                end
            end

            for iBS = 1:nBS
                for iLink = 1:L
                    % Downlink
                    if obj.simulation.simulateDownlink
                        if ~isempty(Links{BS{iBS}.ID, iLink}) && strcmp(Links{BS{iBS}.ID, iLink}.Type, 'Primary') && ~Links{BS{iBS}.ID, iLink}.IsMUSTFarUE
                            Links{BS{iBS}.ID, iLink}.schedule = linkScheduleDL{iBS} == iLink;
                            Links{BS{iBS}.ID, iLink}.firstSubcarrier = find( mean(linkScheduleDL{iBS} == iLink,2), 1, 'first') -1;
                            Links{BS{iBS}.ID, iLink}.scheduledSubcarriers = sum(mean(linkScheduleDL{iBS} == iLink,2));                         
                            if Links{BS{iBS}.ID, iLink}.IsMUSTNearUE
                                FarUEID = Links{BS{iBS}.ID, iLink}.MUSTFarUE;
                                Links{BS{iBS}.ID, FarUEID}.schedule = Links{BS{iBS}.ID, iLink}.schedule;
                                Links{BS{iBS}.ID, FarUEID}.firstSubcarrier = Links{BS{iBS}.ID, iLink}.firstSubcarrier;
                                Links{BS{iBS}.ID, FarUEID}.scheduledSubcarriers = Links{BS{iBS}.ID, iLink}.scheduledSubcarriers;                              
                            end
                        end
                    end
                    % Uplink
                    if obj.simulation.simulateUplink
                        if ~isempty(Links{iLink, BS{iBS}.ID}) && strcmp(Links{iLink, BS{iBS}.ID}.Type, 'Primary')
                            Links{iLink, BS{iBS}.ID}.schedule = linkScheduleUL{iBS} == iLink;
                            Links{iLink, BS{iBS}.ID}.firstSubcarrier = find( mean(linkScheduleUL{iBS} == iLink,2), 1, 'first') -1;
                            Links{iLink, BS{iBS}.ID}.scheduledSubcarriers = sum(mean(linkScheduleUL{iBS} == iLink,2));
                        end
                    end
                end
            end
            
            % assign parameters to desired links
            for linkX = 1:L
                for linkY = 1:L
                    if ~isempty(Links{linkX, linkY})
                        
                        % common objects and properties for primary and interfering links
                        
                        % get base station ID of current link
                        switch Links{linkX, linkY}.direction
                            case 'Downlink'
                                BSiD = Links{linkX, linkY}.Transmitter;
                                %update per BS nonlinearity parameters
                                Links{linkX, linkY}.Nonlinearity = obj.simulation.downlinkNonlinearity;
                                Links{linkX, linkY}.amplifierOBO = obj.simulation.amplifierOBO(BSiD);
                                Links{linkX, linkY}.smoothnessFactor = obj.simulation.smoothnessFactor(BSiD);

                            case 'Uplink'
                                BSiD = Links{linkX, linkY}.Receiver;
                            otherwise
                                error('Link type unkown!');
                        end
                        
                        % set objects and parameters for primary links
                        if strcmp(Links{linkX, linkY}.Type, 'Primary')

                            % get MCS
                            mcs = obj.modulation.mcsValues(Links{linkX, linkY}.ID);

                            % pathloss
                            Links{linkX, linkY}.pathloss            = obj.simulation.pathloss(Links{linkX, linkY}.ID);  
                            
                            % channel coding
                            switch Links{linkX, linkY}.direction
                                case 'Downlink'
                                    maxNCodewords = length(obj.layerMapping.table.Downlink{end});
                                case 'Uplink'
                                    maxNCodewords = length(obj.layerMapping.table.Uplink{end});
                                otherwise 
                                    error('Layer mapping table for this connection unknown');
                            end
                               
                            for iCodeword = 1:maxNCodewords
                                Links{linkX, linkY}.ChannelCoder{iCodeword}        = Coding.ChannelCoding( obj.coding.code{Links{linkX, linkY}.ID}, ... % Coding scheme
                                                                                            obj.coding.decoding{Links{linkX, linkY}.ID}, ...            % Decoding algorithim
                                                                                            mcs.codingRateTimes1024/1024, ...                           % Code rate
                                                                                            obj.coding.decodingIterations(Links{linkX, linkY}.ID) ...   % Decoding iterations
                                                                                            );
                            end
                            
                            % modulation
                            switch obj.modulation.waveform{BSiD}
                                case{'OFDM'}
                                    Links{linkX, linkY}.Modulator   = Modulation.Modulator( 'QAM', mcs.modulationOrder, 'OFDM',...
                                                                                            {'none'},...                                    % data symbol spreading type
                                                                                            obj.simulation.noisePowerEstimation,...         % noise power estimation
                                                                                            Links{linkX, linkY}.scheduledSubcarriers,...    % Number subcarriers
                                                                                            obj.modulation.MCSymbols(BSiD),...              % Number OFDM Symbols
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...      % Subcarrier spacing (Hz)
                                                                                            obj.modulation.samplingRate,...                 % Sampling rate (Samples/s)
                                                                                            Links{linkX, linkY}.firstSubcarrier*...         % Intermediate frequency first subcarrier (Hz)
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...
                                                                                            false,...                                       % Transmit real valued signal
                                                                                            obj.modulation.nGuardSymbols(BSiD)/...          % Cyclic prefix length (s)
                                                                                            (obj.modulation.subcarrierSpacing(BSiD)*...
                                                                                            obj.modulation.MCSymbols(BSiD)), ...  
                                                                                            0 ...                                           % Zero guard length (s)
                                                                                            );
                                case{'f-OFDM'}
                                    Links{linkX, linkY}.Modulator   = Modulation.Modulator( 'QAM', mcs.modulationOrder, 'f-OFDM', ...
                                                                                            {'none'},...                                    % data symbol spreading type
                                                                                            obj.simulation.noisePowerEstimation,...         % noise power estimation
                                                                                            Links{linkX, linkY}.scheduledSubcarriers,...    % Number subcarriers
                                                                                            obj.modulation.MCSymbols(BSiD),...              % Number OFDM Symbols
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...      % Subcarrier spacing (Hz)
                                                                                            obj.modulation.samplingRate,...                 % Sampling rate (Samples/s)
                                                                                            Links{linkX, linkY}.firstSubcarrier*...         % Intermediate frequency first subcarrier (Hz)
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...
                                                                                            false,...                                       % Transmit real valued signal
                                                                                            0,...                                           % Cyclic prefix length (s)
                                                                                            0, ...                                          % Zero guard length (s)
                                                                                            1.5*obj.modulation.nGuardSymbols(BSiD)/...      % Length of the transmit filter (s)
                                                                                            ((obj.modulation.nSymbolsTotal(BSiD)- ...
                                                                                            obj.modulation.nGuardSymbols(BSiD))*...
                                                                                            obj.modulation.subcarrierSpacing(BSiD)),...
                                                                                            1.5*obj.modulation.nGuardSymbols(BSiD)/...      % Length of the receive filter (s)
                                                                                            (obj.modulation.MCSymbols(BSiD)*...
                                                                                            obj.modulation.subcarrierSpacing(BSiD)),...
                                                                                            obj.modulation.nGuardSymbols(BSiD)/...          % Length of the additional cyclic prefix (s)
                                                                                            (obj.modulation.MCSymbols(BSiD)*...
                                                                                            obj.modulation.subcarrierSpacing(BSiD))...
                                                                                            );
                                case{'WOLA'}
                                    Links{linkX, linkY}.Modulator   = Modulation.Modulator( 'QAM', mcs.modulationOrder, 'WOLA', ...
                                                                                            {'none'},...                                    % data symbol spreading type
                                                                                            obj.simulation.noisePowerEstimation,...         % noise power estimation
                                                                                            Links{linkX, linkY}.scheduledSubcarriers,...    % Number subcarriers
                                                                                            obj.modulation.MCSymbols(BSiD),...              % Number OFDM Symbols
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...      % Subcarrier spacing (Hz)
                                                                                            obj.modulation.samplingRate,...                 % Sampling rate (Samples/s)
                                                                                            Links{linkX, linkY}.firstSubcarrier*...         % Intermediate frequency first subcarrier (Hz)
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...
                                                                                            false,...                                       % Transmit real valued signal
                                                                                            0,...                                           % Cyclic prefix length (s)
                                                                                            0,...                                           % Zero guard length (s)
                                                                                            0.5*obj.modulation.nGuardSymbols(BSiD)/...      % Length of the window overlapping (s) at the transmitter
                                                                                            (obj.modulation.MCSymbols(BSiD)*...
                                                                                            obj.modulation.subcarrierSpacing(BSiD)),...
                                                                                            0.5*obj.modulation.nGuardSymbols(BSiD)/...      % Length of the window overlapping (s) at the receiver
                                                                                            (obj.modulation.MCSymbols(BSiD)*...
                                                                                            obj.modulation.subcarrierSpacing(BSiD))...
                                                                                            );
                                case{'FBMC'}
                                   Links{linkX, linkY}.Modulator    = Modulation.Modulator( 'PAM', sqrt(mcs.modulationOrder), 'FBMC',...
                                                                                            {'none'},...                                    % data symbol spreading type
                                                                                            obj.simulation.noisePowerEstimation,...         % noise power estimation
                                                                                            Links{linkX, linkY}.scheduledSubcarriers,...    % Number subcarriers
                                                                                            obj.modulation.MCSymbols(BSiD),...              % Number OFDM Symbols
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...      % Subcarrier spacing (Hz)
                                                                                            obj.modulation.samplingRate,...                 % Sampling rate (Samples/s)
                                                                                            Links{linkX, linkY}.firstSubcarrier*...         % Intermediate frequency first subcarrier (Hz)
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...
                                                                                            false,...                                       % Transmit real valued signal
                                                                                            obj.modulation.prototypeFilter, ...             % prototype filter
                                                                                            4, ...                                          % overlapping factor
                                                                                            0, ...
                                                                                            true...
                                                                                            );
                                case{'UFMC'}
                                    Links{linkX, linkY}.Modulator   = Modulation.Modulator( 'QAM', mcs.modulationOrder, 'UFMC', ...
                                                                                            {'none'},...                                    % data symbol spreading type
                                                                                            obj.simulation.noisePowerEstimation,...         % noise power estimation
                                                                                            Links{linkX, linkY}.scheduledSubcarriers,...    % Number subcarriers
                                                                                            obj.modulation.MCSymbols(BSiD),...              % Number UFMC Symbols
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...      % Subcarrier spacing (Hz)
                                                                                            obj.modulation.samplingRate,...                 % Sampling rate (Samples/s)
                                                                                            Links{linkX, linkY}.firstSubcarrier*...         % Intermediate frequency first subcarrier (Hz)
                                                                                            obj.modulation.subcarrierSpacing(BSiD),...
                                                                                            false,...                                       % Transmit real valued signal
                                                                                            obj.modulation.nSubcarriersPerSubband(BSiD),... % number of subcarriers per subband
                                                                                            round(obj.modulation.nGuardSymbols(BSiD)/...    % Zero prefix length in samples
                                                                                            (obj.modulation.subcarrierSpacing(BSiD)*...
                                                                                            obj.modulation.MCSymbols(BSiD))*obj.modulation.samplingRate),...
                                                                                            0,...                                           % no overhead
                                                                                            1.5*obj.modulation.nGuardSymbols(BSiD)/...      % Length of the filter (s)
                                                                                            (obj.modulation.MCSymbols(BSiD)*...
                                                                                            obj.modulation.subcarrierSpacing(BSiD))...
                                                                                            );
                                otherwise
                                    error('Waveform unkown!');
                            end
                            
                            % assign spatial correlation according to TS36.101, Annex B
                            if strcmp(obj.channel.spatialCorrelation,'none')
                                % no spatial correlation
                                correlationCoefficientTX = 0;
                                correlationCoefficientRX = 0;
                            else
                                switch Links{linkX, linkY}.direction
                                    case 'Downlink'
                                        % for downlink the BS is the TX
                                        correlationCoefficientTX = obj.channel.spatialCorrBS;
                                        correlationCoefficientRX = obj.channel.spatialCorrUE;
                                    case 'Uplink'
                                        % for uplink the UE is the TX
                                        correlationCoefficientTX = obj.channel.spatialCorrUE;
                                        correlationCoefficientRX = obj.channel.spatialCorrBS;
                                    otherwise
                                        error('Link type unkown!');
                                end
                            end
                            
                            % generate channel object for primary links
                            Links{linkX, linkY}.Channel             = Channel.FastFading(   obj.modulation.samplingRate,...
                                                                                            obj.channel.powerDelayProfile,...
                                                                                            Links{linkX, linkY}.Modulator.WaveformObject.Nr.SamplesTotal,...
                                                                                            Links{linkX, linkY}.Velocity / obj.constants.SPEED_OF_LIGHT * obj.simulation.centerFrequency,...
                                                                                            obj.channel.dopplerModel,...
                                                                                            obj.channel.nPaths,...
                                                                                            obj.channel.timeCorrelation,...
                                                                                            obj.channel.spatialCorrelation,...
                                                                                            correlationCoefficientTX,...
                                                                                            correlationCoefficientRX,...
                                                                                            Links{linkX, linkY}.nTxAntennas,...
                                                                                            Links{linkX, linkY}.nRxAntennas,...
                                                                                            false,...
                                                                                            obj.channel.K,...
                                                                                            obj.channel.delta...
                                                                                            );

                            % MIMO
                            Links{linkX, linkY}.MIMO = MIMO.MimoProcessing(obj.modulation.transmissionMode);
                            switch obj.modulation.transmissionMode
                                case {'TxD','CLSM','OLSM'}
                                        codebook                                = Parameters.getCodebook(Links{linkX, linkY}.nTxAntennas,obj.modulation.transmissionMode,obj.modulation.delayDiversity);
                                        Links{linkX, linkY}.MIMO.codebook       = codebook;
                                        if ~isempty(Links{linkX, linkY}.Feedback)
                                            Links{linkX, linkY}.Feedback.codebook = codebook;
                                        end
                                case 'custom'
                                    if ~obj.feedback.enable || ~obj.feedback.pmi 
                                        if size(obj.modulation.precodingMatrix{Links{linkX, linkY}.ID},1) ~= Links{linkX, linkY}.nTxAntennas || size(obj.modulation.precodingMatrix{Links{linkX, linkY}.ID},2)~= obj.modulation.nStreams(Links{linkX, linkY}.ID) 
                                            error('Precoding matrix for link %d does not fit number of antennas and number of streams',Links{linkX, linkY}.ID);
                                        end
                                    end
                                        Links{linkX, linkY}.MIMO.setPrecodingMatrix(obj.modulation.precodingMatrix{Links{linkX, linkY}.ID});
                                        codebook                                = Parameters.getCodebook(Links{linkX, linkY}.nTxAntennas,'CLSM',0);
                                        if ~isempty(Links{linkX, linkY}.Feedback)
                                            Links{linkX, linkY}.Feedback.codebook = codebook;
                                        end
                            end
                            
                            Links{linkX, linkY}.layerMappingTable = obj.layerMapping.table;

                        elseif strcmp(Links{linkX, linkY}.Type, 'Interference')
                            Links{linkX, linkY}.ChannelCoder = [];
                            Links{linkX, linkY}.Modulator = [];
                            Links{linkX, linkY}.schedule = [];
                        end % if link type
                    end % if isempty link
                end % for linkY
            end % for linkX
            
            % assign channel objects to interfering links 
            for linkX = 1:L
                for linkY = 1:L
                    if ~isempty(Links{linkX, linkY}) && strcmp(Links{linkX, linkY}.Type, 'Interference')

                        % find primary BS and get total number of samples
                        switch Links{linkX, linkY}.direction
                            case 'Downlink'
                                for iLink = 1:L
                                    if ~isempty(Links{linkX, iLink}) && strcmp(Links{linkX, iLink}.Type, 'Primary')
                                        nSamplesTotalTmp = Links{linkX, iLink}.Modulator.WaveformObject.Nr.SamplesTotal;
                                    end
                                end
                            case 'Uplink'
                                for iLink = 1:L
                                    if ~isempty(Links{iLink, linkY}) && strcmp(Links{iLink, linkY}.Type, 'Primary')
                                        nSamplesTotalTmp = Links{iLink, linkY}.Modulator.WaveformObject.Nr.SamplesTotal;
                                    end
                                end
                            otherwise
                                error('Link type unkown!');
                        end

                        % generate channel object for interfering links
                        Links{linkX, linkY}.Channel             = Channel.FastFading(   obj.modulation.samplingRate,...
                                                                                        obj.channel.powerDelayProfile,...
                                                                                        nSamplesTotalTmp,...
                                                                                        Links{linkX, linkY}.Velocity / obj.constants.SPEED_OF_LIGHT * obj.simulation.centerFrequency,...
                                                                                        obj.channel.dopplerModel,...
                                                                                        obj.channel.nPaths,...
                                                                                        obj.channel.timeCorrelation,...
                                                                                        obj.channel.spatialCorrelation,...
                                                                                        correlationCoefficientTX,...
                                                                                        correlationCoefficientRX,...
                                                                                        Links{linkX, linkY}.nTxAntennas,...
                                                                                        Links{linkX, linkY}.nRxAntennas,...
                                                                                        false,...
                                                                                        obj.channel.K,...
                                                                                        obj.channel.delta...
                                                                                        );
                    end
                end
            end
            
            % check the generated links
            obj.checkLinks(Links);
            
        end %generateLinks
        
        function checkLinks(obj, Links)
            % this method checks if parameters of the generated links are
            % consistent and correct
            
            L = length(Links);
            
            % check if all modulation schemes (in all cells) approximately
            % fit to the same frame duration
            
            linkCounter = 1;
            for linkX = 1:L
                for linkY = 1:L
                    if ~isempty(Links{linkX, linkY}) && ~strcmp(Links{linkX, linkY}.Type, 'Interference')
                        % for all primary links
                        
                        % save frame duration and symbol time spacing
                        frameDuration(linkCounter,1)    = Links{linkX, linkY}.Modulator.WaveformObject.PHY.TimeSpacing * Links{linkX, linkY}.Modulator.WaveformObject.Nr.MCSymbols;
                        timeSpacing(linkCounter,1)      = Links{linkX, linkY}.Modulator.WaveformObject.PHY.TimeSpacing;
                        
                        linkCounter = linkCounter + 1;
                    end
                end
            end

            
            % check if the difference from maximum to minimum frame
            % duration of all cells is small compared to the symbol time
            if abs(max(frameDuration) - min(frameDuration)) > 0.5 * min(timeSpacing)
                error('Frame durations of cells do not match!');
            end
            
            % save mean frame duration
            obj.simulation.averageFrameDuration = mean(frameDuration);
            
        end %checkLinks
    
    end %methods
    
    
    
end %classdef

