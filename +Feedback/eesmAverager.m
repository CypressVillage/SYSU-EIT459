classdef eesmAverager < handle
    % Implements an EESM averager. Given a set of predefined \beta values
    % and the Modulation and Coding Scheme (MCS) used, returns an effective
    % SINR value.
    % (c) Josep Colom Ikuno, INTHFT, 2008
    % adopted for Link Level by Stefan Schwarz 2010

   properties
       % Place where we will store the beta value for each possible MCS
       beta_tables
       % MCSs corresponding to the previous beta values
       MCS_values
       % Indicate which MCSs are valid
       isvalid
   end

   methods
       function obj = eesmAverager(beta_values,MCS_values)
           % Input checking
           if length(beta_values) ~= length(MCS_values)
               error('The vector containg the beta values and the vector specifying the corresponging MCSs are of different length');
           end
           if min(MCS_values)<0
               error('The minimum MCS cannot be lower than 0');
           end
           
           % Translate the MCSs so we can put them in a lookup table.
           for i_=1:length(beta_values)
               MCS_value_idx = MCS_values(i_) ;
               obj.MCS_values(MCS_value_idx)  = MCS_values(i_);
               obj.beta_tables(MCS_value_idx) = beta_values(i_);
               obj.isvalid(MCS_value_idx)     = true;
           end
       end
       
       % Average the given SINR vector. varargin contains the following:
       %   - varargin{1} = MCS -> values in the range 0:15
       % INPUT VALUES IN LINEAR
       % OUTPUT VALUE IN dB
       % Allows for calculations with multiple beta values
       function effective_SINR = average(obj,SINR_vector,varargin)
           % Ensure matrix dimensions
           SINR_vector = SINR_vector(:);
           MCS_idx = obj.MCS_values;
           effective_SINR = zeros(length(MCS_idx),1);
           if min(MCS_idx)<1
               error('CQI cannot be lower than 0');
           elseif max(MCS_idx)>length(obj.isvalid)+1
               error('CQI cannot be higher than %d',length(obj.isvalid));
           elseif sum(~obj.isvalid(MCS_idx))~=0
               error('SINR averaging not defined for all of the CQIs');
           end
           betas = obj.beta_tables(MCS_idx);
           % Needed due to some strange handling of vectors from Matlab:
           % SINR_vector(ja) results in a column vector, while betas(jb) a
           % row.
           if length(betas)>1 && length(SINR_vector)>1
               % Change to allow to calculate EESMs for multiple betas in a
               % single pass (needed for the scheduler)
               na = length(SINR_vector);
               nb = length(betas);
               [ja,jb] = meshgrid(1:na,1:nb);
               effective_SINR(:) = -(betas.').*log(mean(exp(-SINR_vector(ja)./betas(jb)),2));
           elseif length(SINR_vector)==1
               effective_SINR(:) = -betas.*log(exp(-SINR_vector./betas));
           else
               effective_SINR(:) = -betas.*log(mean(exp(-SINR_vector./betas)));
           end
     
           effective_SINR = 10*log10(effective_SINR);
           effective_SINR(effective_SINR == -Inf) = -450;
       end
   end
end 
