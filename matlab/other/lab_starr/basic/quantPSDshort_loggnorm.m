function [norm_rest_beta norm_active_beta subfreq] = quantPSDshort_loggnorm(rest,active,FREQ_QPSD,order,freq,filename)
% quantPSD performs quantitative analysis on rest/active structure
% arrays with the given frequency ranges. 

% Created by SAS 4/29/10 for ALC's project
% Normalizes M1 power in each of low and high beta bands to high gamma power
% (76-100Hz), as a normalization factor for "background" activity.
%
% called by: Ecog_GroupPSDabsS1norm.m
%
% output: one .mat file with two 3D matrices and # ecog contact that is
% closest to M1.

% SAS 4/29/10: section below commented out because allfreq analysis is not
% relevant to log power gamma normalization analysis
% ---------------------------------------------
% % initialize and populate the 2x3x4 matrix 'allfreq'
% % allfreq is a 2x3x4 3-D matrix in which the four 2x3 arrays correspond to the
% % 4 data recording channels (premotor ecog, M1 ecog, M1-S1 ecog, and stn
% % lfp). **2/19/09 - focusing on M1,S1 now, neglecting premotor for now -ALC
% % The 2 rows contain:
% %   1st row     -   resting state data
% %   2nd row     -   active movement data
% % The 3 columns contain:
% %   1st col     -   frequency at which max power occurs
% %   2nd col     -   max log power value
% %   3rd col     -   total log power across all 5 frequency bands
% 
% allfreq = zeros(2,3,3);
% 
% for i=1:3
%     if isnan(order(i))
%         allfreq(:,:,i) = NaN;
%         continue;
%     end
%     %     rest_data = rest((order(i)),:); % from matrix-based code
%     %     active_data = active((order(i)),:);
%     rest_data = rest.contact_pair(order(i)).log_mean_PSD;
%     active_data = active.contact_pair(order(i)).log_mean_PSD;
%     % Need to limit rest and active data to those freq range of interest (typically 4Hz-100Hz)
%     %the _tmp variables are to limit the frequency spectrum for allfreq max
%     %values without affecting the later subfreq code, which also calls the
%     %rest, active, and freq variables
%     
%     %to make more flexible, will use linear indexing to not require graphs to
%     %have 5 frequency bins
%     %     rest_data_tmp = rest_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
%     %     active_data_tmp = active_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
%     %     freq_tmp = freq(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
%     
%     %Need to change the _tmp such that the max freq is based on the freq bins
%     %and excludes, for example, 55-65Hz.
%     rest_data_tmp = rest_data(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));
%     active_data_tmp = active_data(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));
%     freq_tmp = freq(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));
%     
%     [c1 i1] = max(rest_data_tmp);
%     allfreq(1,1,i)=freq_tmp(i1);
%     allfreq(1,2,i)=c1;
%     [c2 i2] = max(active_data_tmp);
%     allfreq(2,1,i)=freq_tmp(i2);
%     allfreq(2,2,i)=c2;
%     array1 = [];
%     array2 = [];
%     for j = 1:size(FREQ_QPSD,1)
%         tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
%         tmp2 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
%         array1 = [array1 tmp1];  %#ok<AGROW>
%         array2 = [array2 tmp2]; %#ok<AGROW>
%     end
%     allfreq(1,3,i)=sum(array1);
%     allfreq(2,3,i)=sum(array2);
% end
% ---------------------------------------------





% initialize and populate the 5x5x3 matrix 'subfreq'
% The four 5x5 arrays of subfreq correspond to the 3 data channels
% (M1,S1,STN LFP).  Each of the 5 rows correspond to the 5
% frequency bands defined by variable FREQ_QPSD.
% The 5 columns contain:
%   1st col     -   mean log power in given frequency band at rest
%   2nd col     -   col1 normalized by mean log power in high gamma band at rest
%   3rd col     -   mean log power in given frequency band during movement
%   4th col     -   col3 normalized by mean log power in high gamma at active
%   5th col     -   normalized mean log power during movement - mean log power at rest (col 4 - col 2)
for i = 1:size(FREQ_QPSD,1)
    subfreq = zeros(i,5,3);
end

norm_rest_beta=[];
norm_active_beta=[];

for i=1:3
    
    if isnan(order(i))
        subfreq(:,:,i)=NaN;
        continue;
    end

    rest_data = rest.contact_pair(order(i)).log_mean_PSD;
    active_data = active.contact_pair(order(i)).log_mean_PSD;
    
    % pull out log power of relevant frequency bands
    rest_lo_beta = rest_data(freq>FREQ_QPSD(2,1) & freq<FREQ_QPSD(2,2));
    rest_hi_beta = rest_data(freq>FREQ_QPSD(3,1) & freq<FREQ_QPSD(3,2));
    rest_beta = [rest_lo_beta rest_hi_beta];
    rest_gamma = rest_data(freq>FREQ_QPSD(5,1) & freq<FREQ_QPSD(5,2));
    
    active_lo_beta = active_data(freq>FREQ_QPSD(2,1) & freq<FREQ_QPSD(2,2));
    active_hi_beta = active_data(freq>FREQ_QPSD(3,1) & freq<FREQ_QPSD(3,2));
    active_beta = [active_lo_beta active_hi_beta];
    active_gamma = active_data(freq>FREQ_QPSD(5,1) & freq<FREQ_QPSD(5,2));
    
    % normalize beta power by gamma power
    norm_rest_beta = [norm_rest_beta mean(rest_beta)-mean(rest_gamma)];
    norm_active_beta = [norm_active_beta mean(active_beta)-mean(active_gamma)];
    
    for j=1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp2 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        subfreq(j,1,i) = mean(tmp1);
        subfreq(j,2,i) = mean(tmp1)-mean(rest_gamma);
        subfreq(j,3,i) = mean(tmp2);
        subfreq(j,4,i) = mean(tmp2)-mean(active_gamma);
        subfreq(j,5,i) = subfreq(j,4,i)-subfreq(j,2,i);
    end
    
end
% save('allfreq','subfreq');
% %% Save and write quantPSD data as excel spreadsheet
% [data4export] = exportdata(allfreq, subfreq, fn);
% cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.
return;