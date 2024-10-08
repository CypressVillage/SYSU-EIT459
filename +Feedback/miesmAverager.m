classdef miesmAverager < handle
    % Implements an MIESM averager. Given a set of predefined \beta values
    % and the Modulation and Coding Scheme (MCS) used, returns an effective
    % SINR value.
    % (c) Stefan Schwarz 2010


   properties
       % Place where we will store the beta value for each possible MCS
       beta_tables
       beta_tables_dB
       % MCSs corresponding to the previous beta values
       MCS_values
       % Indicate which MCSs are valid
       isvalid
       % BICM data
       MI_data
   end

   methods
       function obj = miesmAverager(beta_values,MCS_values)
           % Input checking
           if length(beta_values) ~= length(MCS_values)
               error('The vector containg the beta values and the vector specifying the corresponging MCSs are of different length');
           end
           if min(MCS_values)<0
               error('The minimum MCS cannot be lower than 0');
           end
           
           % Translate the MCSs so we can put them in a lookup table.
           for i_=1:length(beta_values)
               MCS_value_idx = MCS_values(i_);
               obj.MCS_values(MCS_value_idx)  = MCS_values(i_);
               obj.beta_tables(MCS_value_idx) = beta_values(i_);
               obj.isvalid(MCS_value_idx)     = true;
           end
           MI_data_loader = Feedback.MIdata_loader;
           obj.MI_data = MI_data_loader.MIdata_load;
           
           % Calculation in the log domain directly (faster)
           for i_=1:length(obj.MI_data)
               obj.MI_data(i_).SNR_dB = 10*log10(obj.MI_data(i_).SNR);
               
               % We assume that the SNRs are equidistant in dBs!! If not... the behavior can be unexpected.
               obj.MI_data(i_).SNR_dB_step = mean(obj.MI_data(i_).SNR_dB(2:end)-obj.MI_data(i_).SNR_dB(1:(end-1)));
           end
           obj.beta_tables_dB = 10*log10(obj.beta_tables);
       end
       
       % Average the given SINR vector. varargin contains the following:
       %   - varargin{1} = MCS -> values in the range 0:15
       %   - varargin{2} = symbol alphabet order corresponding to the MCS values
       % INPUT VALUES IN LINEAR
       % OUTPUT VALUE IN dB
       % Allows for calculations with multiple beta values
       function effective_SINR_dB = average(obj,SINR_vector,varargin)
           % Ensure matrix dimensions
           SINR_vector    = SINR_vector(:);
           MCS_idx        = obj.MCS_values;
           alphabets      = log2(varargin{1})/2;
           effective_SINR = zeros(length(MCS_idx),1);
           
           
           
           % Optional mode in which the input is directly in dB
           if length(varargin)<3
               input_in_dB = false;
           else
               input_in_dB = varargin{3};
           end
           
           if min(MCS_idx)<0
               error('CQI cannot be lower than 0');
           elseif max(MCS_idx)>length(obj.isvalid)
               error('CQI cannot be higher than %d',length(obj.isvalid));
           elseif sum(~obj.isvalid(MCS_idx))~=0
               error('SINR averaging not defined for all of the CQIs');
           end
           
           % Check for upper ceiling
           if ~input_in_dB
               larger = (SINR_vector > obj.MI_data(3).SNR(end));
%                larger = (SINR_vector > obj.MI_data(3).SNR(19));
           else               
               larger = (SINR_vector > obj.MI_data(3).SNR_dB(end));
%                larger = (SINR_vector > obj.MI_data(3).SNR_dB(19));
           end
           if sum(abs(diff(SINR_vector))) > 0
               if sum(larger) == length(SINR_vector) % in this case all SINRs are larger than the saturation value of BICMC and therefore we simply average them...
                   if ~input_in_dB
                       effective_SINR = mean(SINR_vector);
                   else
                       effective_SINR = mean(10.^(SINR_vector/10));
                   end
               elseif length(SINR_vector) > 1
                   for cqi_i = MCS_idx
                       cqi_idx = cqi_i;
                       cqi_MI_data = obj.MI_data(alphabets(cqi_idx));
                       
                       % if cqi_i
                       if ~input_in_dB
                           SINR_vector_temp = SINR_vector/obj.beta_tables(cqi_idx);
                       else
                           SINR_vector_temp = SINR_vector-obj.beta_tables_dB(cqi_idx);
                       end
                       % SINR_vector_temp(SINR_vector_temp > obj.MI_data(alphabets(cqi_i+1)).max_SNR) = obj.MI_data(alphabets(cqi_i+1)).max_SNR;
                       % SINR_vector_temp(SINR_vector_temp < obj.MI_data(alphabets(cqi_i+1)).min_SNR) = obj.MI_data(alphabets(cqi_i+1)).min_SNR;
                       if ~input_in_dB
                           ind = zeros(1,length(SINR_vector_temp));
                           for SNR_i = 1:length(SINR_vector_temp)
                               [~,ind(SNR_i)] = min(abs(SINR_vector_temp(SNR_i)-cqi_MI_data.SNR));
                           end
                       else
                           dB_step     = cqi_MI_data.SNR_dB_step;
                           init_SNR_dB = cqi_MI_data.SNR_dB(1);
                           max_ind     = length(cqi_MI_data.SNR_dB);
                           ind         = round((SINR_vector_temp-init_SNR_dB)/dB_step+1);
                           ind         = max(ind,1);
                           ind         = min(ind,max_ind);
                       end
                       
                       %                    I_avg = mean(min(cqi_MI_data.BICM(ind),5.55));
                       %                    I_avg = mean(cqi_MI_data.BICM(ind));
                       I_avg = mean(cqi_MI_data.BICM(ind))-std(cqi_MI_data.BICM(ind))*0;%0.15; % for large SNR variations, the mean overestimates the performance of LTE --> we slightly bias the averaging towards smaller values
                       [~,ind2] = min(abs(I_avg-cqi_MI_data.BICM));
                       if ~input_in_dB
                           effective_SINR(cqi_idx) = obj.beta_tables(cqi_idx)*cqi_MI_data.SNR(ind2);
                       else
                           effective_SINR(cqi_idx) = obj.beta_tables_dB(cqi_idx)+cqi_MI_data.SNR_dB(ind2);
                       end
                       % I_avg = mean(interp1(obj.MI_data(alphabets(cqi_i)).SNR,obj.MI_data(alphabets(cqi_i)).BICM,SINR_vector_temp,'linear','extrap'));
                       % effective_SINR(cqi_i+1) = obj.beta_tables(cqi_i)*interp1(obj.MI_data(alphabets(cqi_i)).BICM(1:obj.MI_data(alphabets(cqi_i)).sat_SNR),obj.MI_data(alphabets(cqi_i)).SNR(1:obj.MI_data(alphabets(cqi_i)).sat_SNR),I_avg,'linear','extrap');
                       % else
                       %     effective_SINR(cqi_i+1) = 10^-10;
                       % end
                   end
               else
                   effective_SINR = SINR_vector;
               end
           else
               effective_SINR = SINR_vector(1);
           end
           
           % Final conversion (if necessary)
           if ~input_in_dB
               effective_SINR_dB = 10*log10(effective_SINR);
           else
               effective_SINR_dB = effective_SINR;
           end
       end
   end
end
