function [subfreq] = quantDesync(rest,active,FREQ_QPSD,order,freq,filename)
% quantPSD performs quantitative analysis on rest/active structure
% arrays with the given frequency ranges.  allows user to define ecog
% contact closest to M1 **UPDATED 2/19/09 to work with ecog_lfp data -ALC




% initialize and populate the 5x5x4 matrix 'subfreq'
% The four 5x5 arrays of subfreq correspond to the 4 data channels
% (pre-motor,M1,M1-S1,STN LFP).  Each of the 5 rows correspond to the 5
% frequency bands defined by variable FREQ_QPSD.
% The 5 columns contain:
%   1st col     -   total power in given frequency band at rest
%   2nd col     -   power in given frequency band, divided by total power
%                   in all 5 freq bands, at rest
%   3rd col     -   total power in given frequency band during movement
%   4th col     -   power in given frequency band, divided by total power
%                   in all 5 freq bands, during movement
%   5th col     -   power during movement divided by power at rest (col 3
%                   divided by col 1)
for i = 1:size(FREQ_QPSD,1)
    subfreq = zeros(i,5,3);
end

for i=1:3
    if isnan(order(i))
    subfreq(:,:,i)=NaN;
    continue;
    end
    rest_data = rest.contact_pair(order(i)).mean_PSD;
    active_data = active.contact_pair(order(i)).mean_PSD;
    for j=1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp2 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        subfreq(j,1,i) = sum(tmp1);
        subfreq(j,2,i) = sum(tmp1);
        subfreq(j,3,i) = sum(tmp2);
        subfreq(j,4,i) = sum(tmp2);
        subfreq(j,5,i) = sum(tmp2)/sum(tmp1);
    end
end
% save('allfreq','subfreq');
% %% Save and write quantPSD data as excel spreadsheet
% [data4export] = exportdata(allfreq, subfreq, fn);
% cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.
return;