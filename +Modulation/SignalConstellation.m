classdef SignalConstellation < handle 
    % =====================================================================        
    % This MATLAB class represents a QAM or PAM signal constellation. The 
    % parameters are initialized by the class contructor. 
    % Afterwards we can transform a bit stream into symbols by 
    % ".Bit2Symbol()" and back from symbol to bit by ".Symbol2Bit".
    % Furthermore LLR calculation for an AWGN channal and MIMO is also
    % included. 
    % =====================================================================    
    % Ronald Nissel, rnissel@nt.tuwien.ac.at
    % (c) 2017 by Institute of Telecommunications, TU Wien
    % www.tc.tuwien.ac.at
    % =====================================================================    

  properties (SetAccess = private)
       Method
       ModulationOrder
       BitMapping
       SymbolMapping
       Implementation
       text
  end
  
  % LTE-A-Pro Mapping w/ 3GPP MUST parameters
  properties
        BPSK_I = sqrt(1/2)*[1,-1].';
        BPSK_Q = sqrt(1/2)*[1,-1].';
        QPSK_I = sqrt(1/2)*[1,1,-1,-1].';
        QPSK_Q = sqrt(1/2)*[1,-1,1,-1].';
        QAM16_I = sqrt(1/10)*[1,1,3,3,1,1,3,3,-1,-1,-3,-3,-1,-1,-3,-3].';
        QAM16_Q = sqrt(1/10)*[1,3,1,3,-1,-3,-1,-3,1,3,1,3,-1,-3,-1,-3].';
        QAM64_I = sqrt(1/42)*[3,3,1,1,3,3,1,1,5,5,7,7,5,5,7,7,3,3,1,1,3,3,1,1,5,5,7,7,5,5,7,7,-3,-3,-1,-1,-3,-3,-1,-1,-5,-5,-7,-7,-5,-5,-7,-7,-3,-3,-1,-1,-3,-3,-1,-1,-5,-5,-7,-7,-5,-5,-7,-7].';
        QAM64_Q = sqrt(1/42)*[3,1,3,1,5,7,5,7,3,1,3,1,5,7,5,7,-3,-1,-3,-1,-5,-7,-5,-7,-3,-1,-3,-1,-5,-7,-5,-7,3,1,3,1,5,7,5,7,3,1,3,1,5,7,5,7,-3,-1,-3,-1,-5,-7,-5,-7,-3,-1,-3,-1,-5,-7,-5,-7].';
        QAM256_I = sqrt(1/170)*[5,5,7,7,5,5,7,7,3,3,1,1,3,3,1,1,5,5,7,7,5,5,7,7,3,3,1,1,3,3,1,1,11,11,9,9,11,11,9,9,13,13,15,15,13,13,15,15,11,11,9,9,11,11,9,9,13,13,15,15,13,13,15,15,5,5,7,7,5,5,7,7,3,3,1,1,3,3,1,1,5,5,7,7,5,5,7,7,3,3,1,1,3,3,1,1,11,11,9,9,11,11,9,9,13,13,15,15,13,13,15,15,11,11,9,9,11,11,9,9,13,13,15,15,13,13,15,15,-5,-5,-7,-7,-5,-5,-7,-7,-3,-3,-1,-1,-3,-3,-1,-1,-5,-5,-7,-7,-5,-5,-7,-7,-3,-3,-1,-1,-3,-3,-1,-1,-11,-11,-9,-9,-11,-11,-9,-9,-13,-13,-15,-15,-13,-13,-15,-15,-11,-11,-9,-9,-11,-11,-9,-9,-13,-13,-15,-15,-13,-13,-15,-15,-5,-5,-7,-7,-5,-5,-7,-7,-3,-3,-1,-1,-3,-3,-1,-1,-5,-5,-7,-7,-5,-5,-7,-7,-3,-3,-1,-1,-3,-3,-1,-1,-11,-11,-9,-9,-11,-11,-9,-9,-13,-13,-15,-15,-13,-13,-15,-15,-11,-11,-9,-9,-11,-11,-9,-9,-13,-13,-15,-15,-13,-13,-15,-15].';
        QAM256_Q = sqrt(1/170)*[5,7,5,7,3,1,3,1,5,7,5,7,3,1,3,1,11,9,11,9,13,15,13,15,11,9,11,9,13,15,13,15,5,7,5,7,3,1,3,1,5,7,5,7,3,1,3,1,11,9,11,9,13,15,13,15,11,9,11,9,13,15,13,15,-5,-7,-5,-7,-3,-1,-3,-1,-5,-7,-5,-7,-3,-1,-3,-1,-11,-9,-11,-9,-13,-15,-13,-15,-11,-9,-11,-9,-13,-15,-13,-15,-5,-7,-5,-7,-3,-1,-3,-1,-5,-7,-5,-7,-3,-1,-3,-1,-11,-9,-11,-9,-13,-15,-13,-15,-11,-9,-11,-9,-13,-15,-13,-15,5,7,5,7,3,1,3,1,5,7,5,7,3,1,3,1,11,9,11,9,13,15,13,15,11,9,11,9,13,15,13,15,5,7,5,7,3,1,3,1,5,7,5,7,3,1,3,1,11,9,11,9,13,15,13,15,11,9,11,9,13,15,13,15,-5,-7,-5,-7,-3,-1,-3,-1,-5,-7,-5,-7,-3,-1,-3,-1,-11,-9,-11,-9,-13,-15,-13,-15,-11,-9,-11,-9,-13,-15,-13,-15,-5,-7,-5,-7,-3,-1,-3,-1,-5,-7,-5,-7,-3,-1,-3,-1,-11,-9,-11,-9,-13,-15,-13,-15,-11,-9,-11,-9,-13,-15,-13,-15].';
        
        % 3GPP MUST composite constellation parameters
        c_QPSK = [sqrt(1/5),2/sqrt(29),7*sqrt(1/578)];
        d_QPSK = [sqrt(2),5/(2*sqrt(2)),23/(7*sqrt(2))];
        c_QAM16 = [sqrt(5/21),3*sqrt(5/334),sqrt(5/69)];
        d_QAM16 = [2*sqrt(2/5),17/(3*sqrt(10)),4*sqrt(2/5)];
        c_QAM64 = [sqrt(21/85),sqrt(7/34),sqrt(7/55)];
        d_QAM64 = [4*sqrt(2/21),3*sqrt(3/14),2*sqrt(6/7)];
        MUSTNearUEScale
        MUSTIdx = 0;
  end
  
  
  
  methods
  	function obj = SignalConstellation(ModulationOrder,Method, MUSTIdx)
        % Generates a QAM or PAM signal constellation object with 
        % corresponding Gray-coded bit mapping. The first argument 
        % represents the modulation order and the second argument the
        % method, either 'QAM' or 'PAM'. For example (4,'QAM'), (256,'QAM')
        % or (2,'PAM'), (16,'PAM')
        obj.MUSTIdx = MUSTIdx;
        obj.ModulationOrder = ModulationOrder;
        obj.Method          = Method;
        
        % Gray coded bitmapping
        if strcmp(obj.Method,'QAM')
            switch ModulationOrder
                case 2
                    I = obj.BPSK_I;
                    Q = obj.BPSK_Q;
                case 4
                    I = obj.QPSK_I;
                    Q = obj.QPSK_Q;
                    if obj.MUSTIdx ~= 0
                        c = obj.c_QPSK(obj.MUSTIdx);
                        d = obj.d_QPSK(obj.MUSTIdx);
                    end
                case 16
                    I = obj.QAM16_I;
                    Q = obj.QAM16_Q;
                    if obj.MUSTIdx ~= 0
                        c = obj.c_QAM16(obj.MUSTIdx);
                        d = obj.d_QAM16(obj.MUSTIdx);     
                    end
                case 64
                    I = obj.QAM64_I;
                    Q = obj.QAM64_Q;
                    if obj.MUSTIdx ~= 0
                        c = obj.c_QAM64(obj.MUSTIdx);
                        d = obj.d_QAM64(obj.MUSTIdx);
                    end
                case 256
                    I = obj.QAM256_I;
                    Q = obj.QAM256_Q;
            end
            if obj.MUSTIdx == 0
                obj.SymbolMapping = I + 1i*Q;
                obj.BitMapping = de2bi(0:obj.ModulationOrder-1);
            else
                obj.ModulationOrder = obj.ModulationOrder *4;
                It = c*(I-d);
                Qt = c*(Q-d);
                % Inner mapping is NearUE, outter mapping is FarUE (restricted to legacy LTE QPSK mapping).
                obj.SymbolMapping = [-It - 1i*Qt; -It + 1i*Qt; It + -1i*Qt; It + 1i*Qt];
                obj.BitMapping = de2bi(0:obj.ModulationOrder-1);
                obj.MUSTNearUEScale = c;
            end          
        elseif strcmp(obj.Method,'PAM')
            obj.BitMapping = [ones(obj.ModulationOrder/2,1);zeros(obj.ModulationOrder/2,1)];
            for i_temp = 2:log2(obj.ModulationOrder)
                BinaryTemp = obj.BitMapping(1:2:end,i_temp-1);
                obj.BitMapping(:,i_temp) = [BinaryTemp;BinaryTemp(end:-1:1)];
            end
            obj.SymbolMapping = (2*(1:obj.ModulationOrder)-obj.ModulationOrder-1).';
            obj.SymbolMapping = obj.SymbolMapping/sqrt(mean(abs(obj.SymbolMapping).^2));
            % Determine the underlying symbol alphabet and the corresponding
            % bit mapping
            [~,SortOrder] = sort(bi2de(obj.BitMapping),'ascend');
            obj.SymbolMapping = obj.SymbolMapping(SortOrder);
            obj.BitMapping = obj.BitMapping(SortOrder,:);
        else
           error('Signal constellation method must be QAM or PAM!');
        end

        % For the LLR detection, we determine all data symbols which have a
        % bit value of one (zero) at a certain position

        obj.Implementation.DataSymbolsBitvalueOne   = repmat(obj.SymbolMapping,1,log2(obj.ModulationOrder));
        obj.Implementation.DataSymbolsBitvalueOne   = reshape(obj.Implementation.DataSymbolsBitvalueOne(logical(obj.BitMapping)),obj.ModulationOrder/2,[]);
        obj.Implementation.DataSymbolsBitvalueZero  = repmat(obj.SymbolMapping,1,log2(obj.ModulationOrder));
        obj.Implementation.DataSymbolsBitvalueZero  = reshape(obj.Implementation.DataSymbolsBitvalueZero(not(logical(obj.BitMapping))),obj.ModulationOrder/2,[]);       
    end   
    
    function DataSymbols = Bit2Symbol(obj,PrimaryUEStream, MUSTUEStream)
        % Maps a bit stream to the correpsonding symbol alphabet
        tmpSize = size(PrimaryUEStream);
        if obj.MUSTIdx == 0
            DataSymbols = obj.SymbolMapping( bi2de(reshape(PrimaryUEStream(:),log2(obj.ModulationOrder),[])')+1 );
            DataSymbols = reshape( DataSymbols, tmpSize(1)/log2(obj.ModulationOrder), tmpSize(2) );
        else           
            bits = [reshape(PrimaryUEStream, [], log2(obj.ModulationOrder/4)) reshape(MUSTUEStream', 2, [])'];
            DataSymbols = obj.SymbolMapping(bi2de(bits) + 1);
        end
    end

    function EstimatedBitStream = Symbol2Bit(obj,EstimatedDataSymbols)
        % Maps symbols (nearest neighbor detection) to the corresponding
        % bit stream
        EstimatedDataSymbols = EstimatedDataSymbols(:);
        
        [~,b] = min(abs((repmat(EstimatedDataSymbols,1,obj.ModulationOrder)-repmat((obj.SymbolMapping).',size(EstimatedDataSymbols,1),1)).'));
        EstimatedBitStream = obj.BitMapping(b(:),1:log2(obj.ModulationOrder)).';      
        EstimatedBitStream = EstimatedBitStream(:);
    end
    function QuantizedDataSymbols = SymbolQuantization(obj,EstimatedDataSymbols)
        % Performs quantization of the received symbols, that is nearest
        % neighbor detection
        EstimatedDataSymbols = EstimatedDataSymbols(:);

        [~,b] = min(abs((repmat(EstimatedDataSymbols,1,obj.ModulationOrder)-repmat((obj.SymbolMapping).',size(EstimatedDataSymbols,1),1)).'));
        QuantizedDataSymbols = obj.SymbolMapping(b(:),:).';
        QuantizedDataSymbols = QuantizedDataSymbols(:);
    end    
    
    function LLR = LLR_AWGN(obj,y,Pn)
        % Calculates the LLR for an AWGN channel, that is, y=x+n with y
        % denoting the received data symbol, x the transmitted data symbol
        % and n the Gaussian distributed noise with power Pn
    
        if numel(Pn)>1
            PnRepeated = reshape(repmat(Pn.',log2(obj.ModulationOrder)*obj.ModulationOrder/2,1),obj.ModulationOrder/2,[]);
        else 
            PnRepeated = Pn;
        end       
        ReceivedDataSymbolsRepeated     = reshape(repmat(y.',log2(obj.ModulationOrder)*obj.ModulationOrder/2,1),obj.ModulationOrder/2,[]);
        DataSymbolsBitvalueOneRepeated  = repmat(obj.Implementation.DataSymbolsBitvalueOne,1,length(y));
        DataSymbolsBitvalueZeroRepeated = repmat(obj.Implementation.DataSymbolsBitvalueZero,1,length(y));
 
        LLR =   log(sum(exp(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueOneRepeated).^2./PnRepeated),1)./...
                    sum(exp(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueZeroRepeated).^2./PnRepeated),1)).';
        
        if sum(isnan(LLR)) % too large numbers in the exponential --> apply max-log approximation
            LLR = (max(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueOneRepeated).^2./PnRepeated,[],1) - max(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueZeroRepeated).^2./PnRepeated,[],1)).';
        end

        LLR(LLR==Inf)=10^10;
        LLR(LLR==-Inf)=-10^10;
    end
    
    function LLR = LLR_MIMO_ML(obj,y,H,Rn, MIMOmethod, Precoder, PrecodingMatrix)
        % This method is a straighhfoward (high compuational complexity)
        % implementation of a MIMO maximum likelihood LLR calculation (per
        % bit). The model assumes y=Hx+n, whereas "y" represents the received
        % symbols of size y(NrRX,NrPos), "H" the channel matrix of size 
        % H(NrRX,NrTx,NrPos), "Rn" the correlation matrix of the Gaussian 
        % distributed noise with size Rn(NrRX,NrRX,NrPos), and 
        % "Precoder" the precoding matrix of size Precoder(NrTx,NrStreams)

        if not(exist('PrecodingMatrix','var'))
            PrecodingMatrix = eye(size(H,2));
        end
  
        NrTxStreams = size(PrecodingMatrix,2);
        NrTimeFrequencyPositions = size(H,3);
                     
        IndexRef = repmat(1:obj.ModulationOrder,obj.ModulationOrder^(NrTxStreams-1),1);
        IndexRef = IndexRef(:);
        IndexDataSymbolAntennaCombination = zeros(length(IndexRef),NrTxStreams);
        for i=0:NrTxStreams-1
            IndexDataSymbolAntennaCombination(:,i+1) = repmat(IndexRef(1:obj.ModulationOrder^i:end),obj.ModulationOrder^i,1);
        end
        IndexDataSymbolAntennaCombination = IndexDataSymbolAntennaCombination(1:end/2,:)';
                          
        for i_BitPos = 1:size(obj.Implementation.DataSymbolsBitvalueOne,2)
            % Matlab changes dimensions!!!!!
            if NrTxStreams==2
                x_One_Atom = [obj.Implementation.DataSymbolsBitvalueOne(IndexDataSymbolAntennaCombination(1,:),i_BitPos).';obj.SymbolMapping(IndexDataSymbolAntennaCombination(2,:)).'];
                x_Zero_Atom = [obj.Implementation.DataSymbolsBitvalueZero(IndexDataSymbolAntennaCombination(1,:),i_BitPos).';obj.SymbolMapping(IndexDataSymbolAntennaCombination(2,:)).'];
            else
                x_One_Atom = [obj.Implementation.DataSymbolsBitvalueOne(IndexDataSymbolAntennaCombination(1,:),i_BitPos).';obj.SymbolMapping(IndexDataSymbolAntennaCombination(2:end,:))];
                x_Zero_Atom = [obj.Implementation.DataSymbolsBitvalueZero(IndexDataSymbolAntennaCombination(1,:),i_BitPos).';obj.SymbolMapping(IndexDataSymbolAntennaCombination(2:end,:))];            
            end 
            for i_TX_Antenna = 0:NrTxStreams-1                
                x_One(:,:,i_TX_Antenna+1,i_BitPos) = circshift(x_One_Atom,[i_TX_Antenna,0]);
                x_Zero(:,:,i_TX_Antenna+1,i_BitPos) = circshift(x_Zero_Atom,[i_TX_Antenna,0]);  
            end            
        end
        
        [~,b,c,d] = size(x_One);
        a = size(y,1);
        LLR = zeros(c,d,NrTimeFrequencyPositions);
        
        if strcmp(MIMOmethod,'OLSM')
            NrTxStreams = size(Precoder.U,1);
            l = 1:NrTimeFrequencyPositions;
            if size(Precoder.W,1)==4
                k = mod(floor((l-1)/NrTxStreams),4)+1;
            else
                k=ones(1,NrTimeFrequencyPositions);
            end
                
            p = mod(l-1,NrTxStreams)+1;
        end
        
        for i_TimeFrequencyPos = 1:NrTimeFrequencyPositions
            
            Rn_temp = Rn(:,:,i_TimeFrequencyPos);
            InvRn_temp = Rn_temp^-1;
            [S,V] = svd(InvRn_temp);            
            C = (S*sqrt(V))';
            
            if strcmp(MIMOmethod,'OLSM')
               PrecodingMatrix = Precoder.W{k(i_TimeFrequencyPos)}*(Prcoder.D).^(p(i_TimeFrequencyPos))*Precoder.U;
            end
           
            y_temp = C*y(:,i_TimeFrequencyPos);
            H_temp = C*H(:,:,i_TimeFrequencyPos)*PrecodingMatrix;
                        
            y_temp_repeated = repmat(y_temp,1,b,c,d);            
            Hx_One_temp = reshape(H_temp*x_One(:,:),a,b,c,d);
            Hx_Zero_temp = reshape(H_temp*x_Zero(:,:),a,b,c,d);
               
            LLR(:,:,i_TimeFrequencyPos)=squeeze(log(sum(exp(sum(-abs(y_temp_repeated-Hx_One_temp).^2,1)),2)./sum(exp(sum(-abs(y_temp_repeated-Hx_Zero_temp).^2,1)),2)));  
        end
        LLR = LLR(:,:).';
        LLR(LLR==Inf)=10000;
        LLR(LLR==-Inf)=-10000;   
    end  
    
    function [LLR,x_est,NoiseScaling] = LLR_MIMO_ZF(obj,y,H,Pn, MIMOmethod, Precoder, PrecodingMatrix, despreading, nSubcarriers)
        % This method calculates the LLR for a MIMO system after a ZF
        % equalization. To keep the compuational complexity low, the cross
        % correlation after the ZF equalizer is neglected. 
        % The model assumes y=Hx+n, whereas "y" represents the received
        % symbols of size y(NrRX,NrPos), "H" the channel matrix of size 
        % H(NrRX,NrTx,NrPos), "Pn" the power of the Gaussian noise n 
        % and "Precoder" the precoding matrix of size Precoder(NrTx,NrStreams) 
        
        if not(exist('PrecodingMatrix','var'))
            PrecodingMatrix = eye(size(H,2));
        end  

        N = size(y,2);
        NrTxStreams = size(PrecodingMatrix,2);
        NrRxAntennas = size(y,1);
        
        if strcmp(MIMOmethod,'OLSM')
            NrTxStreams = size(Precoder.W{1},2);
            l = 1:N;
            if size(Precoder.W,1)==4
                k = mod(floor((l-1)/NrTxStreams),4)+1;
            else
                k=ones(1,N);
            end
            p = mod(l-1,NrTxStreams)+1;
        elseif strcmp(MIMOmethod,'TxD')
            [H, y] = obj.TxD_decoding(H, y);
            PrecodingMatrix = 1;
            NrTxStreams = size(H,2);
            N = size(y,2);
        end
        
%         N = size(H,3);
        x_est = zeros(NrTxStreams,N);
        NoiseScaling = zeros(NrTxStreams,N);        
        for i_n = 1:N
           if strcmp(MIMOmethod,'OLSM')
               PrecodingMatrix = Precoder.W{k(i_n)}*(Precoder.D).^(p(i_n))*Precoder.U;
           end

           H_temp = H(:,:,i_n)*PrecodingMatrix;
           
           if NrRxAntennas<=NrTxStreams
                ZF_Equalizer = H_temp'/(H_temp*H_temp');
           else
                ZF_Equalizer = (H_temp'*H_temp)^-1*H_temp';
           end  
%            ZF_Equalizer = pinv(H_temp);
           x_est(:,i_n) = (ZF_Equalizer*y(:,i_n));           
           NoiseScaling(:,i_n) = Pn*sum(ZF_Equalizer.*conj(ZF_Equalizer),2);
        end
        
        
        if strcmp(MIMOmethod,'TxD')
            x_est = reshape(x_est,[],1);
            x_est(2:2:end) = conj(x_est(2:2:end));
        else
            x_est = x_est.';
        end
        
        NrTxStreams = size(x_est,2);
        % symbol de-spreading
        for iStream = 1:NrTxStreams
            x_est(:,iStream) = despreading(x_est(:,iStream),nSubcarriers);
        end
        
        NoiseScaling = NoiseScaling.';
        LLR = reshape(obj.LLR_AWGN(reshape(x_est,[],1), reshape(NoiseScaling,[],1)),[],NrTxStreams); 
    end 
    
    
    function [LLR,x_est,NoiseScaling,UnbiasedScaling] = LLR_MIMO_MMSE(obj,y,H,Pn, MIMOmethod, Precoder, PrecodingMatrix, despreading, nSubcarriers)
        % This method calculates the LLR for a MIMO system after an 
        % unbiased MMSE equalization. To keep the compuational complexity
        % low, the cross correlation after the MMSE equalizer is neglected. 
        % The model assumes y=Hx+n, whereas "y" represents the received
        % symbols of size y(NrRX,NrPos), "H" the channel matrix of size 
        % H(NrRX,NrTx,NrPos), "Pn" the power of the Gaussian noise n 
        % and "Precoder" the precoding matrix of size Precoder(NrTx,NrStreams) 
        
        if not(exist('PrecodingMatrix','var'))
            PrecodingMatrix = eye(size(H,2));
        end
        
        N = size(y,2);
        NrTxStreams = size(PrecodingMatrix,2);
        NrRXAntennas = size(H,1); 
        
        if strcmp(MIMOmethod,'OLSM')
            NrTxStreams = size(Precoder.U,1);
            l = 1:N;
            if size(Precoder.W,1)==4
                k = mod(floor((l-1)/NrTxStreams),4)+1;
            else
                k=ones(1,N);
            end
            
            p = mod(l-1,NrTxStreams)+1;
        elseif strcmp(MIMOmethod,'TxD')
            [H, y] = obj.TxD_decoding(H, y);
            PrecodingMatrix = 1;
            NrTxStreams = size(H,2);
            [NrRXAntennas, N] = size(y);
        end
        
        x_est = zeros(NrTxStreams,N);
        NoiseScaling = zeros(NrTxStreams,N); 
        UnbiasedScaling = zeros(NrTxStreams,N); 
        for i_n = 1:N
            
           if strcmp(MIMOmethod,'OLSM')
                PrecodingMatrix = Precoder.W{k(i_n)}*(Prcoder.D).^(p(i_n))*Precoder.U;
           end
           
           H_temp = H(:,:,i_n)*PrecodingMatrix;
           MMSE_Equalizer = H_temp'/(H_temp*H_temp'+Pn*eye(NrRXAntennas));  
           x_est(:,i_n) = (MMSE_Equalizer*y(:,i_n));
           
           Temp = MMSE_Equalizer*H_temp;
           NoiseScaling(:,i_n) = Pn*sum(MMSE_Equalizer.*conj(MMSE_Equalizer),2)+sum(abs(Temp-diag(diag(Temp))).^2,2);
           UnbiasedScaling(:,i_n) = abs(sum(MMSE_Equalizer.*H_temp.',2));
        end

        % symbol de-spreading        
        if strcmp(MIMOmethod,'TxD')
            x_est = reshape(x_est,[],1);
            x_est(2:2:end) = conj(x_est(2:2:end));
            NoiseScaling = reshape(NoiseScaling,[],1).';
            UnbiasedScaling = reshape(UnbiasedScaling,[],1).';
        else
            x_est = x_est.';
        end
        NrTxStreams = size(x_est,2);
        
        for iStream = 1:NrTxStreams
            x_est(:,iStream) = despreading(x_est(:,iStream),nSubcarriers);
        end
        
        NoiseScaling = NoiseScaling.';
        UnbiasedScaling = UnbiasedScaling.';
        LLR = reshape(obj.LLR_AWGN(reshape(x_est./UnbiasedScaling,[],1), reshape(NoiseScaling./UnbiasedScaling.^2,[],1)),[],NrTxStreams); 
    end   
    
    function LLR = LLR_MIMO_Sphere(obj,y,H,Pn, MIMOmethod, Precoder, PrecodingMatrix)
        % This method calculates the LLR for a MIMO system by employing
        % a SphereDecoder. The MATLAB Communications Toolbox is REQUIRED!
        % The model assumes y=Hx+n, whereas "y" represents the received
        % symbols of size y(NrRX,NrPos), "H" the channel matrix of size 
        % H(NrRX,NrTx,NrPos), "Pn" the power of the Gaussian noise n 
        % and "Precoder" the precoding matrix of size Precoder(NrTx,NrStreams)
        NrTxStreams = size(PrecodingMatrix,2);

        if strcmp(MIMOmethod,'OLSM')
            NrTxStreams = size(Precoder.U,1);
            l = 1:size(H,3);
            if size(Precoder.W,1)==4
                k = mod(floor((l-1)/NrTxStreams),4)+1;
            else
                k=ones(1,size(H,3));
            end
            p = mod(l-1,NrTxStreams)+1;
        end
        
        if exist('PrecodingMatrix','var')
            Heff = nan(size(H,1),NrTxStreams,size(H,3));
            for i_H = 1:size(H,3)
                if strcmp(MIMOmethod,'OLSM')
                    PrecodingMatrix = Precoder.W{k(i_H)}*(Prcoder.D).^(p(i_H))*Precoder.U;
                end
                Heff(:,:,i_H)=H(:,:,i_H)*PrecodingMatrix;
            end
        else
            Heff = H;
        end
        
        Heff_permute = permute(Heff,[3 2 1]);
        CommToolboxSphereDecoder = comm.SphereDecoder('Constellation', obj.SymbolMapping,'BitTable', double(obj.BitMapping), 'DecisionType', 'Soft');
        LLR = step(CommToolboxSphereDecoder,y.',Heff_permute)*(1/Pn);       
    end
    
    
    function [ H_combined, r_combined ] = TxD_decoding(obj, H, r)
            [nRx, nTx, NrPos] = size(H);
            H = permute(H, [3,1,2]);
            r = permute(r, [2,1]);
            switch nTx

                case 2
                    switch nRx
                        case 2
                            r_combined = zeros( 4, NrPos/2);
                            H_combined = zeros( 4, nTx, NrPos/2);
                             for j = 1:(NrPos/2)
                                    i = 2*(j-1) + 1;
                                    r_combined(:,j) = [r(i,1);r(i,2);conj(r(i+1,1));conj(r(i+1,2))];
                                    H_combined(:,:,j) = [H(i,1,1),-H(i,1,2);H(i,2,1),-H(i,2,2);conj(H(i+1,1,2)),conj(H(i+1,1,1));conj(H(i+1,2,2)),conj(H(i+1,2,1))];
                             end
                        case 1
                            r_combined = zeros( 2, NrPos/2);
                            H_combined = zeros( 2, nTx, NrPos/2);
                             for j = 1:(NrPos/2)
                                    i = 2*(j-1) + 1;
                                    r_combined(:,j) = [r(i);conj(r(i+1))];
                                    H_combined(:,:,j) = [H(i,1,1),-H(i,1,2);conj(H(i+1,1,2)),conj(H(i+1,1,1))]; 
                             end
                        otherwise
                        error('The mode transmit diversity with %d transmit antennas is not implemented for %d receive antennas.', nTx, nRx);
                    end

                case 4
                    switch nRx
                        case 4
                            r_combined = zeros( 16, NrPos/4);
                            H_combined = zeros( 16, nTx, NrPos/4);
                            for j = 1:(NrPos/4)
                                i = 4*(j-1) + 1;
                                r_combined(:,j) = [r(i,1);r(i,2);conj(r(i+1,1));conj(r(i+1,2));
                                                r(i+2,1);r(i+2,2);conj(r(i+3,1));conj(r(i+3,2));
                                                r(i,3);r(i,4);conj(r(i+1,3));conj(r(i+1,4));
                                                r(i+2,3);r(i+2,4);conj(r(i+3,3));conj(r(i+3,4))];

                                H_combined(:,:,j) = [H(i,1,1),-H(i,1,3),0,0;H(i,2,1),-H(i,2,3),0,0;conj(H(i+1,1,3)),conj(H(i+1,1,1)),0,0;conj(H(i+1,2,3)),conj(H(i+1,2,1)),0,0;...  % generate equivalent channel matrix
                                    0,0,H(i+2,1,2),-H(i+2,1,4);0,0,H(i+2,2,2),-H(i+2,2,4);0,0,conj(H(i+3,1,4)),conj(H(i+3,1,2));0,0,conj(H(i+3,2,4)),conj(H(i+3,2,2));...
                                    H(i,3,1),-H(i,3,3),0,0;H(i,4,1),-H(i,4,3),0,0;conj(H(i+1,3,3)),conj(H(i+1,3,1)),0,0;conj(H(i+1,4,3)),conj(H(i+1,4,1)),0,0;...  
                                    0,0,H(i+2,3,2),-H(i+2,3,4);0,0,H(i+2,4,2),-H(i+2,4,4);0,0,conj(H(i+3,3,4)),conj(H(i+3,3,2));0,0,conj(H(i+3,4,4)),conj(H(i+3,4,2))];
                            end
                        case 2
                            r_combined = zeros( 8, NrPos/4);
                            H_combined = zeros( 8, nTx, NrPos/4);
                            for j = 1:(NrPos/4)
                                i = 4*(j-1) + 1;
                                r_combined(:,j) = [r(i,1);r(i,2);conj(r(i+1,1));conj(r(i+1,2));
                                                r(i+2,1);r(i+2,2);conj(r(i+3,1));conj(r(i+3,2))];

                                H_combined(:,:,j) = [H(i,1,1),-H(i,1,3),0,0;H(i,2,1),-H(i,2,3),0,0;conj(H(i+1,1,3)),conj(H(i+1,1,1)),0,0;conj(H(i+1,2,3)),conj(H(i+1,2,1)),0,0;...  % generate equivalent channel matrix
                                    0,0,H(i+2,1,2),-H(i+2,1,4);0,0,H(i+2,2,2),-H(i+2,2,4);0,0,conj(H(i+3,1,4)),conj(H(i+3,1,2));0,0,conj(H(i+3,2,4)),conj(H(i+3,2,2))];
                            end

                        case 1
                            r_combined = zeros( 4, NrPos/4);
                            H_combined = zeros( 4, nTx, NrPos/4);
                            for j = 1:(NrPos/4)
                                i = 4*(j-1) + 1;
                                r_combined(:,j) = [r(i,1);conj(r(i+1,1));r(i+2,1);conj(r(i+3,1))];

                                H_combined(:,:,j) = [H(i,1,1),-H(i,1,3),0,0;conj(H(i+1,1,3)),conj(H(i+1,1,1)),0,0;
                                    0,0,H(i+2,1,2),-H(i+2,1,4);0,0,conj(H(i+3,1,4)),conj(H(i+3,1,2))];
                            end
                        otherwise
                        error('The mode transmit diversity with %d transmit antennas is not implemented for %d receive antennas.', nTx, nRx);
                    end

                case 8
                    switch nRx
                        case 2 
                            r_combined = zeros( 16, NrPos/8);
                            H_combined = zeros( 16, nTx, NrPos/8);
                            for j = 1:(NrPos/8)
                                i = 8*(j-1) + 1;
                                r_combined(:,j) = [r(i,1);r(i,2);conj(r(i+1,1));conj(r(i+1,2));
                                                r(i+2,1);r(i+2,2);conj(r(i+3,1));conj(r(i+3,2));...
                                                r(i+4,1);r(i+4,2);conj(r(i+5,1));conj(r(i+5,2));...
                                                r(i+6,1);r(i+6,2);conj(r(i+7,1));conj(r(i+7,2))];
                               H_combined(:,:,j) = [H(i,1,1),-H(i,1,5),0,0,0,0,0,0;H(i,2,1),-H(i,2,5),0,0,0,0,0,0;conj(H(i+1,1,5)),conj(H(i+1,1,1)),0,0,0,0,0,0;conj(H(i+1,2,5)),conj(H(i+1,2,1)),0,0,0,0,0,0;...
                                                    0,0,H(i+2,1,2),-H(i+2,1,6),0,0,0,0;0,0,H(i+2,2,2),-H(i+2,2,6),0,0,0,0;0,0,conj(H(i+3,1,6)),conj(H(i+3,1,2)),0,0,0,0;0,0,conj(H(i+3,2,6)),conj(H(i+3,2,2)),0,0,0,0;...
                                                    0,0,0,0,H(i+4,1,3),-H(i+4,1,7),0,0;0,0,0,0,H(i+4,2,3),-H(i+4,2,7),0,0;0,0,0,0,conj(H(i+5,1,7)),conj(H(i+5,1,3)),0,0;0,0,0,0,conj(H(i+5,2,7)),conj(H(i+5,2,3)),0,0;...
                                                    0,0,0,0,0,0,H(i+6,1,4),-H(i+6,1,8);0,0,0,0,0,0,H(i+6,2,4),-H(i+6,2,8);0,0,0,0,0,0,conj(H(i+7,1,8)),conj(H(i+7,1,4));0,0,0,0,0,0,conj(H(i+7,2,8)),conj(H(i+7,2,4))];            


                            end

                        case 1

                            r_combined = zeros( 8, NrPos/8);
                            H_combined = zeros( 8, nTx, NrPos/8);
                            for j = 1:(NrPos/8)
                                i = 8*(j-1) + 1;
                                r_combined(:,j) = [r(i,1);conj(r(i+1,1));r(i+2,1);conj(r(i+3,1));...
                                                    r(i+4,1);conj(r(i+5,1));r(i+6,1);conj(r(i+7,1))];
                               H_combined(:,:,j) = [H(i,1,1),-H(i,1,5),0,0,0,0,0,0;conj(H(i+1,1,5)),conj(H(i+1,1,1)),0,0,0,0,0,0;...
                                                    0,0,H(i+2,1,2),-H(i+2,1,6),0,0,0,0;0,0,conj(H(i+3,1,6)),conj(H(i+3,1,2)),0,0,0,0;...
                                                    0,0,0,0,H(i+4,1,3),-H(i+4,1,7),0,0;0,0,0,0,conj(H(i+5,1,7)),conj(H(i+5,1,3)),0,0;...
                                                    0,0,0,0,0,0,H(i+6,1,4),-H(i+6,1,8);0,0,0,0,0,0,conj(H(i+7,1,8)),conj(H(i+7,1,4))];
                           end

                        otherwise
                        error('The mode transmit diversity with %d transmit antennas is not implemented for %d receive antennas.', nTx, nRx);
                    end
                otherwise
                    error('The mode transmit diversity is not implemented for %d transmit antennas.',nTx);




            end
      end
  end
end
      