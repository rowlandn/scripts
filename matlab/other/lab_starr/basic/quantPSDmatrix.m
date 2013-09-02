function [allfreq subfreq order] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename)
% quantPSD performs quantitative analysis on rest/active structure
% arrays with the given frequency ranges.  allows user to define ecog
% contact closest to M1 **UPDATED 2/19/09 to work with ecog_lfp data -ALC
%
% output: one .mat file with two 3D matrices and # ecog contact that is
% closest to M1.
 
%UPDATED 2/19/09 to work with ecog_lfp data -ALC
%updated 4/6/09 such that allfreq now excludes data less than 4Hz or 
%greater than 100Hz
%%
% select ecog contact closest to M1 and use that to reassign ecog contact
% data array
ncontact = {'1' '2' '3' '4' '5' '6'};
contactm1 = menu('Select ecog contact closest to M1',ncontact);
% assign contact numbers for each structure relative to m1 contact
m1 = contactm1;
% pre = contactm1+1;
% m1s1 = contactm1-2; % updated 2/13:from M1 contact #, subtract 2 instead of 1 to find S1
s1 = contactm1 - 2;
if contactm1==1
    s1 = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
elseif contactm1==2
    s1 = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
% elseif contactm1==5
%     pre = NaN;
%     menu(['There is no contact pair over premotor.' sprintf('\n')...
%         'Click OK to continue'],'OK');
elseif contactm1==6
    m1 = NaN;
    pre = NaN;
    menu(['There is no contact pair over M1 or premotor.' sprintf('\n')...
        'Click OK to continue'],'OK');
end

% look for LFP channel
% num_contact_pair = length(rest.contact_pair);
[num_row num_col]=size(rest); %Need to keep raw data matrix here b/c has rest and active onsets
num_contact_pair = num_row; %last two columns are rest and active onsets
if num_contact_pair == 6
    lfp = 6; % assume that 6th contact_pair in rest/active structures always contains LFP data
elseif num_contact_pair<6
    lfp = NaN;
end

% order = [pre m1 m1s1 lfp];
order = [m1 s1 lfp];

% initialize and populate the 2x3x4 matrix 'allfreq'
% allfreq is a 2x3x4 3-D matrix in which the four 2x3 arrays correspond to the
% 4 data recording channels (premotor ecog, M1 ecog, M1-S1 ecog, and stn
% lfp). **2/19/09 - focusing on M1,S1 now, neglecting premotor for now -ALC
% The 2 rows contain:
%   1st row     -   resting state data
%   2nd row     -   active movement data
% The 3 columns contain:
%   1st col     -   frequency at which max power occurs
%   2nd col     -   max power value
%   3rd col     -   total power across all 5 frequency bands

allfreq = zeros(2,3,3);
for i=1:3
    if isnan(order(i))
        allfreq(:,:,i) = NaN;
        continue;
        
    end
    rest_data = rest((order(i)),:);%selects appropriate row of data corresponding to contact pair being evaluated
    active_data = active((order(i)),:);%as above
% % Need to limit rest and active data to those freq range of interest (typically 4Hz-100Hz)
%     rest_data = rest_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
%     active_data = active_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
%     freq = freq(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
    rest_data_tmp = rest_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
    active_data_tmp = active_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
    freq_tmp = freq(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
    [c1 i1] = max(rest_data_tmp);%finds max value of rest data and that value's index position
    allfreq(1,1,i)=freq_tmp(i1);%puts freq corresponding to max value into allfreq matrix
    allfreq(1,2,i)=c1;%puts max value into allfreq matrix
    [c2 i2] = max(active_data_tmp);
    allfreq(2,1,i)=freq_tmp(i2);
    allfreq(2,2,i)=c2;
    array1 = [];
    array2 = [];
    for j = 1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp2 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        array1 = [array1 tmp1];  %#ok<AGROW>
        array2 = [array2 tmp2]; %#ok<AGROW>
    end
    allfreq(1,3,i)=sum(array1);
    allfreq(2,3,i)=sum(array2);
end

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
subfreq = zeros(5,5,3);
for i=1:3
    if isnan(order(i))
        subfreq(:,:,i)=NaN;
        continue;
    end
    rest_data = rest((order(i)),:);
    active_data = active((order(i)),:);
    for j=1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp2 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        subfreq(j,1,i) = sum(tmp1);
        subfreq(j,2,i) = sum(tmp1)/allfreq(1,3,i);
        subfreq(j,3,i) = sum(tmp2);
        subfreq(j,4,i) = sum(tmp2)/allfreq(2,3,i);
        subfreq(j,5,i) = sum(tmp2)/sum(tmp1);
    end
end
save([filename '_ecogPSD'],'allfreq','subfreq','order');
return;