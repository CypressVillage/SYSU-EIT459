classdef Feedback < handle
    %% The Feedback Class
    %   Mariam Mussbah
    %   
    %   (c) 2017 Institute of Telecommunications, TU Wien.
    %   www.nt.tuwien.ac.at
    %% The Feedback Properties
    properties
        delay                   % Feedback delay
        Pmi                     % Precoding Matrix Indicator
        Ri                      % Rank Indicator
        Cqi                     % Channel Quality Indicator
        nRxAntennas             % Number of receive antenna
        nTxAntennas             % Number of transmit antenna
        codebook                % Code Book
        Csi                     % Full Channel Matrix
        defaultPrecodingMatrix
    end
    %% 
    methods
        function obj = Feedback(feedback, nTxAntennas, nRxAntennas,...
                                 defaultPrecodingMatrix)
            % This is the class constructor
            obj.nRxAntennas             = nRxAntennas;
            obj.nTxAntennas             = nTxAntennas;
            obj.delay                   = feedback.delay;
            obj.Pmi.enable              = feedback.pmi;
            obj.Ri.enable               = feedback.ri;
            obj.Cqi.enable              = feedback.cqi;
            obj.defaultPrecodingMatrix  = defaultPrecodingMatrix;
            if obj.delay==0
                arraysize=obj.delay+1;
            else
                arraysize=obj.delay;
            end
            obj.Pmi.pmiArray            = ones(arraysize,1);
            obj.Ri.riArray              = ones(arraysize,1);
            obj.Csi                     = cell(1,arraysize);
            
            if feedback.cqi 
                obj.Cqi.cqiArray = cell(1,arraysize);
                obj.Cqi.cqiArray(1,1:feedback.delay+1) = {ones(1,1)};
                switch feedback.averager.Type
                case 'eesm'
                    obj.Cqi.averager = Feedback.eesmAverager(feedback.averager.EESMbetas,...
                                                            feedback.averager.MCS_values);
                case 'miesm'
                    obj.Cqi.averager = Feedback.miesmAverager(feedback.averager.MIESMbetas,...
                                                            feedback.averager.MCS_values);
                otherwise
                    error('Averager type not supported!');
                end
                obj.Cqi.mappingTable = feedback.CqiMappingTable;
            end
            obj.Ri.max = min(nTxAntennas,nRxAntennas);
                        
        end
        

        function updateFeedback(obj, channel, pilotMatrix, noisePower, modulationOrder)
            % updateFeedback this function updates the PMI, RI and CSI.
            
            if ~isempty(channel)
%           -------> for full CSI
%                 obj.Csi(1:end-1)=obj.Csi(2:end);
%                 obj.Csi{end}=channel; 
                % channel on pilot position. used for CLSM
                scheduledChannel = reshape( channel( repmat(logical( sum(pilotMatrix,3)),...
                                             1,1, obj.nRxAntennas,  obj.nTxAntennas) ), [], ...
                                            obj.nRxAntennas, obj.nTxAntennas  );
               % channel used for RI calcultion in case of OLSM
%                scheduledChannel = reshape( channel( repmat(~logical( sum(pilotMatrix,3)),...
%                                              1,1, obj.nRxAntennas,  obj.nTxAntennas) ), [], ...
%                                             obj.nRxAntennas, obj.nTxAntennas  );
                
                if obj.Pmi.enable
                    [pmi,ri] = obj.getPmi(scheduledChannel, noisePower);

                    obj.Pmi.pmiArray(1:end-1) = obj.Pmi.pmiArray(2:end);
                    obj.Pmi.pmiArray(end) = pmi;
                    
                    obj.Ri.riArray(1:end-1) = obj.Ri.riArray(2:end);
                    obj.Ri.riArray(end) = ri; 
                end % update pmi
                if ~obj.Pmi.enable && obj.Ri.enable
                    ri = obj.getRI(scheduledChannel, noisePower);
                    obj.Ri.riArray(1:end-1) = obj.Ri.riArray(2:end);
                    obj.Ri.riArray(end) = ri; 
                    
                end
                if obj.Cqi.enable
                    if obj.Pmi.enable
                        cqi = obj.getCqi(scheduledChannel, obj.codebook.W{obj.Pmi.pmiArray(end)},...
                                       noisePower, modulationOrder);
                    elseif obj.Ri.enable
                        cqi = obj.getCqi(scheduledChannel, [], noisePower,...
                                       modulationOrder);
                    else
                        cqi = obj.getCqi(scheduledChannel, obj.defaultPrecodingMatrix, noisePower,...
                                        modulationOrder);
                    end
                    obj.Cqi.cqiArray(1:end-1) = obj.Cqi.cqiArray(2:end);
                    obj.Cqi.cqiArray{end} = cqi;
                end % update cqi
            end 
        end
        
        function [ precodingMatrix ] = getPrecodingMatrix(obj, index)
            if index == 0
                precodingMatrix=obj.defaultPrecodingMatrix;
            elseif index <= size(obj.codebook.W, 2)
                precodingMatrix=obj.codebook.W{index};
            else
                warning('index exceeds codebook size');
            end
        end
    end % public methods
    methods (Access= private)

        function [maxIndex,Ri] = getPmi(obj, scheduledChannel, noisePower)
            % getPmi calculates the PMI for the given channel matrix
            % 'scheduledChannel' and the noise Power 'noisePower'

            codebookSize = size(obj.codebook.W, 2);            
            maxMutualInformation = 0;
            maxIndex = 1;
            for index = 1:codebookSize
                 precodingMatrix = obj.codebook.W{index};
                 if rank(precodingMatrix) <= obj.Ri.max
                        mutualInformation = 0;
                        for iRE = 1:size(scheduledChannel,1)
                            singleEffectiveChannel = reshape(scheduledChannel(iRE,:,:),...
                                        [obj.nRxAntennas,obj.nTxAntennas])* precodingMatrix;
                            equalizer = pinv(singleEffectiveChannel' * singleEffectiveChannel)...
                                        * singleEffectiveChannel';
                            K = equalizer*singleEffectiveChannel;
                            for iLayer = 1:size(K,1)
                                a = sum(abs(K(iLayer,:)).^2)-abs(K(iLayer,iLayer)); % interstream interference
                                b = noisePower*sum(abs(equalizer(iLayer,:))); % noise enhancement
                                SINR = abs(K(iLayer,iLayer))^2/(a+b);
                                mutualInformation = mutualInformation+log2(1+SINR);
                            end
                        end % iter over all RE
                        if(mutualInformation > maxMutualInformation)
                            maxMutualInformation = mutualInformation;
                            maxIndex = index;
                        end                    
                end
            end %iter over all codebooks
            Ri = rank(obj.codebook.W{maxIndex});
        end %getPMI
        
        
        function [cqi] = getCqi(obj, scheduledChannel, precodingMatrix, noisePower,...
                                modulationOrder)
            % getCqi calculates the CQI for the given channel matrix
            % 'scheduledChannel' and the noise Power 'noisePower'
            N = size(scheduledChannel, 1);
            if isempty(precodingMatrix)
                RI = obj.Ri.riArray(end);
                l = 1:N;
                if size(obj.codebook.W,1)==4
                    k = mod(floor((l-1)/RI),4)+1;
                else
                    k=ones(1,N);
                end
                
                p = mod(l-1,RI)+1;
                OLSM = true;
                
            else
                RI = rank(precodingMatrix);
                OLSM = false;
            end
            sinr = zeros(N, RI);
            for iRE = 1:N
                if OLSM
                   precodingMatrix = obj.codebook.W{k(iRE),RI}*(obj.codebook.D{RI}).^(p(iRE))*obj.codebook.U{RI};
                end
                singleEffectiveChannel = reshape(scheduledChannel(iRE,:,:), [obj.nRxAntennas,...
                                                obj.nTxAntennas])* precodingMatrix;
                equalizer = pinv(singleEffectiveChannel' * singleEffectiveChannel)*...
                                    singleEffectiveChannel';
                K = equalizer*singleEffectiveChannel;
                for l = 1:size(precodingMatrix,2)
                        a = sum(abs(K(l,:)).^2)-abs(K(l,l)); %interstream interference
                        b = noisePower*sum(abs(equalizer(l, :)).^2); %noise
                        sinr(iRE,l) = abs(K(l, l))^2/(a+b);
                end %layers
            end
            cqi = obj.mapSnrToCqi(sinr, modulationOrder); 
        end %getCQI
        
        function cqi = mapSnrToCqi(obj, sinr, modulationOrder)
            % mapSnrToCqi maps the SINR vector of each layer to an
            % equivalent effective SNR and then returns the corresponding CQI
            
            cqi = zeros(obj.Ri.riArray(1), 1);
            if size(sinr, 2)>1
                for l=1:size(sinr, 3)
                    effSNR = obj.Cqi.averager.average(sinr(:,:,l), modulationOrder);
                    temp = zeros(size(obj.Cqi.averager.MCS_values));
                    temp(obj.Cqi.mappingTable(obj.Cqi.averager.MCS_values) <= effSNR) = 1;
                    cqi(l) = find(temp, 1, 'last');
                    if isempty(cqi(l))
                        cqi(l) = 1;
                    end
                end
            else
                effSNR = obj.Cqi.averager.average(sinr, modulationOrder);
                temp = zeros(size(obj.Cqi.averager.MCS_values));
                temp(obj.Cqi.mappingTable(obj.Cqi.averager.MCS_values) <= effSNR) = 1;
                cqi = find(temp,1,'last');
                if isempty(cqi)
                    cqi = 1;
                end
            end
        end % mapSnrToCqi
        
         function maxRI = getRI(obj, scheduledChannel, noisePower)
                    
            maxMutualInformation = 0;
            maxRI = 1;
            l = 1:size(scheduledChannel,1);
            
            for RI = 1:obj.Ri.max
                        if size(obj.codebook.W,1)==4
                            k = mod(floor((l-1)/RI),4)+1;
                        else
                            k=ones(1,size(scheduledChannel,1));
                        end
                            p = mod(l-1,RI)+1;

                        mutualInformation = 0;
                        for n = 1:size(scheduledChannel,2)
                            precodingMatrix = obj.codebook.W{k(n),RI}*(obj.codebook.D{RI}).^(p(n))*obj.codebook.U{RI};
                            singleEffectiveChannel = reshape(scheduledChannel(n,:,:),...
                                        [obj.nRxAntennas,obj.nTxAntennas])* precodingMatrix;
                            equalizer = pinv(singleEffectiveChannel' * singleEffectiveChannel)...
                                        * singleEffectiveChannel';
                            K = equalizer*singleEffectiveChannel;
                            for iLayer = 1:size(K,1)
                                a = sum(abs(K(iLayer,:)).^2)-abs(K(iLayer,iLayer)); % interstream interference
                                b = noisePower*sum(abs(equalizer(iLayer,:))); % noise enhancement
                                SINR = abs(K(iLayer,iLayer))^2/(a+b);
                                mutualInformation = mutualInformation+log2(1+SINR);
                            end
                        end % end n
                        
                        if mutualInformation > maxMutualInformation
                            maxMutualInformation = mutualInformation;
                            maxRI = RI;
                        end                    
                
            end % end iter RI 

         end % end getRI
        
    end %methods
    
    
end

