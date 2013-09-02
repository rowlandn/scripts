function [allfreq subfreq order] = quantPSDrest(psdall,FREQ_QPSD,freq,filename,M1_ch)
% quantPSD performs quantitative analysis on "true rest" data (output = psdall variable
% from ecogPSD_restdata program.  allows user to define ecog
% contact closest to M1
%
% output: one .mat file with two 3D matrices and # ecog contact that is
% closest to M1.

%% define M1, S1, LFP
m1 = M1_ch;
% pre = contactm1+1;
% m1s1 = contactm1-2; % updated 2/13:from M1 contact #, subtract 2 instead of 1 to find S1
s1 = m1 - 2; %updated 6/9/09 for clarity

if  m1==1
    s1 = NaN;
    menu(['There is no contact pair over S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
elseif m1==2
    s1 = NaN;
    menu(['There is no contact pair over S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
% elseif contactm1==5
%     pre = NaN;
%     menu(['There is no contact pair over premotor.' sprintf('\n')...
%         'Click OK to continue'],'OK');
elseif m1==6
    m1 = NaN;
%     pre = NaN;
    menu(['There is no contact pair over M1 or premotor.' sprintf('\n')...
        'Click OK to continue'],'OK');
end

% look for LFP channel
num_contact_pair = size(psdall,2);
if num_contact_pair == 6
    lfp = 6; % assume that 6th contact_pair in rest/active structures always contains LFP data
elseif num_contact_pair<6
    lfp = NaN;
else
    lfp = NaN; % added b/c epilepsy data may have num_contact_pair>6
end

% order = [pre m1 m1s1 lfp];
order = [m1 s1 lfp];

F = freq;
%% initialize and populate the 1x3x3 matrix 'allfreq'
% allfreq is a 1x3x3 3-D matrix in which the 3 arrays correspond to the
% 3 data recording channels ( M1 ecog, S1 ecog, and stn lfp). 
% The row contains data from "true rest" condition:
% The 3 columns contain:
%   1st col     -   frequency at which max power occurs
%   2nd col     -   max power value
%   3rd col     -   total power across all 5 frequency bands

allfreq = zeros(1,3,3);
for i=1:3
    if isnan(order(i))
        allfreq(:,:,i) = NaN;
        continue;
    end
    rest_data = psdall(:,(order(i)));
    F_tmp = F(F>FREQ_QPSD(1,1) & F<FREQ_QPSD(5,2));
    rest_data_tmp = rest_data(F>FREQ_QPSD(1,1) & F<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
    [c1 i1] = max(rest_data_tmp);
    allfreq(1,1,i)=F_tmp(i1);
    allfreq(1,2,i)=c1;
    array1 = [];
    for j = 1:size(FREQ_QPSD,1)
        tmp1 = rest_data(F>FREQ_QPSD(j,1) & F<FREQ_QPSD(j,2));
        array1 = [array1; tmp1];  %#ok<AGROW>
    end
    allfreq(1,3,i)=sum(array1);
%     allfreq(1,3,i) = sum(rest_data_tmp);
end

%% initialize and populate the matrix 'subfreq'
% The 5x2x3 arrays of subfreq correspond to the 4 data channels
% (pre-motor,M1,M1-S1,STN LFP).  Each of the 5 rows correspond to the 5
% frequency bands defined by variable FREQ_QPSD.
% The 5 columns contain:
%   1st col     -   total power in given frequency band
%   2nd col     -   mean power in given frequency band, divided by sum of
%                   mean powers from 5 frequency bands.  Percentage should
%                   add up to 100%.
%   3rd col     -   mean power in given frequency band minus mean power in
%                   high gamma.  Used in gamma-normalized power at rest
subfreq = zeros(5,3,3);
for i=1:3
    if isnan(order(i))
        subfreq(:,:,i)=NaN;
        continue;
    end
    rest_data = psdall(:,(order(i)));

    for j=1:size(FREQ_QPSD,1)
        tmp1 = rest_data(F>FREQ_QPSD(j,1) & F<FREQ_QPSD(j,2));
        subfreq(j,1,i) = mean(tmp1);
    end
    
    %Mean power for each frequency band divided by sum mean power across 5 frequency bands
    for j=1:size(FREQ_QPSD,1)
        subfreq(j,2,i) = subfreq(j,1,i)/sum(subfreq(:,1,i));
        subfreq(j,3,i) = subfreq(j,1,i)-subfreq(5,1,i);
    end
end
% save([filename '_ecogPSD'],'allfreq','subfreq','order');
return;