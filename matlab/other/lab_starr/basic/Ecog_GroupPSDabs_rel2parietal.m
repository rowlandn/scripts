function Ecog_GroupPSDabs_rel2parietal
% Definition: same as Ecog_GroupPSD, except bars represent the actual power
% in each frequency band, rather than the percent power relative to the
% total power across the five frequency bands - ALC 3/9/10

% INPUT = group _ecogPSD files from quantPSD data folder
% specifically, needs mean PSD data

% CALLS function 'quantPSDshort' - this is nearly identical to the
% quantPSD subfunction called by ecogPSD, except that it does not assign
% contact numbers or create the variable 'order' - as it already exists.
% The purpose of calling 'quantPSDshort'is to recalcuate the allfreq and
% subfreq variables in the event that the user wishes to change the
% FREQ_QPSD variable for the purpose of analyzing group data.
% 'quantPSDshort' also does not write data to an excel file.

% OUTPUT: none, creates plot

% 3/9/10: this is identical to Ecog_GroupPSD except instead of graphing the
% % power of each freq band (subfreq(:,2) for rest and subfreq(:,4)) for
% active state, it is looking for "absolute" power (subfreq(:,1) and
% subfreq(:,3) respectively)**ALC
%% define variables for plotting
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   
% FREQ_QPSD = [4 12;...     % delta alpha band
%              13 30;...   % low beta band
%              31 55;...   % high beta band
%              76 100];    % high gamma band     
% FREQBANDS = {num2str(FREQ_QPSD(1,:)), num2str(FREQ_QPSD(2,:)) num2str(FREQ_QPSD(3,:))... 
%     num2str(FREQ_QPSD(4,:)) num2str(FREQ_QPSD(5,:))};
% AREAS = {'pre-motor' 'M1' 'S1' 'STN LFP'};
BRAINAREAS = {'M1' 'S1' 'STN LFP'};
LIMBAREAS = {'HAND' 'ELBOW' 'SHOULDER' 'JAW' 'FOOT' '"ARM"' '"NON-ARM"'};
% YLIM = [0 80];
%% Define frequency bands based on FREQ_QPSD 
%allows for user to create desired number of freq bins
FREQBANDS = [];
for i = 1:size(FREQ_QPSD,1)
    FREQBANDS = [FREQBANDS {num2str(FREQ_QPSD(i,:))}];
end
%% choose channel
k = menu('Select brain area for analysis',BRAINAREAS);
l = menu('Select limb for analysis',LIMBAREAS);

%% grab PD data
% pathnamePD = uigetdir('','Select directory with PD ecogPSD.mat files');
pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numPD = length(PDdir);

PDrest = [];
PDactive=[];
counterPD = 0;

for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    if isnan(order(k));
        continue; % skip files without contact pair over selected area
    end
    
    % Due to brain atrophy, sulci, vessels, etc, some subjects experieinced
% strong ECOG signals, while others had weaker signals. This creates a fair
% amount of variability in the data, given our small group sizes. To
% address this we will calculate the strength of the M1 signal relative to
% the parietal signal (contacts 1v2).
   for j = 1:size(rest.contact_pair,2)
      
    rest.contact_pair(1,j).rel_signal = rest.contact_pair(1,j).mean_PSD ./ ...
        mean(rest.contact_pair(1,1).mean_PSD);
    active.contact_pair(1,j).rel_signal = active.contact_pair(1,j).mean_PSD ...
        ./ mean(active.contact_pair(1,1).mean_PSD);
   end
    
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    [allfreq subfreq] = quantPSDshort_rel2parietal(rest,active,FREQ_QPSD,order,freq,filename);
    
    %now create matrix of percent power in each freq band for each subject
    %to calculate "absolute" power in each band (rather than "% power") use
    %subfreq(:,1,k) and subfreq(:,3,k) for rest and active respectively
    PDrest = [PDrest subfreq(:,1,k)];%#ok<AGROW> %2nd column contains log power at rest
    PDactive = [PDactive subfreq(:,3,k)]; %#ok<AGROW> %4th column contains log power with movement
    counterPD = counterPD+1;
end

Yrest(:,1) = mean(PDrest,2);
% Erest(:,1) = std(PDrest,0,2);
Erest(:,1) = std(PDrest,0,2)/(sqrt(counterPD));
Yactive(:,1) = mean(PDactive,2);
% Eactive(:,1) = std(PDactive,0,2);
Eactive(:,1) = std(PDactive,0,2)/(sqrt(counterPD));

%% Calculate log power relative to each subject's overall signal strength
 


%% grab dystonia data
% pathnameDYT = uigetdir('','Select directory with dystonia ecogPSD.mat files');
dyspath = uigetdir('', 'Select directory that contains Dys _ecgPSD files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
% FileListDYS = dir('*_ecogPSD.mat');

% nfilesDYT = length(FileListDYT);
numDYS = length(DYSdir);

DYSrest=[];
DYSactive=[];
counterDYS=0;
for i=1:numDYS
    filename = DYSdir(i).name;
    load(filename);
    if isnan(order(k))
        continue; %skip files without contact pair over selected area
    end
    
     % Due to brain atrophy, sulci, vessels, etc, some subjects experieinced
% strong ECOG signals, while others had weaker signals. This creates a fair
% amount of variability in the data, given our small group sizes. To
% address this we will calculate the strength of the M1 signal relative to
% the parietal signal (contacts 1v2).
   for j = 1:size(rest.contact_pair,2)
    rest.contact_pair(1,j).rel_signal = rest.contact_pair(1,j).mean_PSD ./ ...
        mean(rest.contact_pair(1,1).mean_PSD);
    active.contact_pair(1,j).rel_signal = active.contact_pair(1,j).mean_PSD ...
        ./ mean(active.contact_pair(1,1).mean_PSD);
   end
    
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    [allfreq subfreq] = quantPSDshort_rel2parietal(rest,active,FREQ_QPSD,order,freq,filename);
    
    
    
    %now create matrix of percent power in each freq band for each subject
    %to calculate "absolute" power in each band (rather than "% power") use
    %subfreq(:,1,k) and subfreq(:,3,k) for rest and active respectively
    DYSrest = [DYSrest subfreq(:,1,k)]; %#ok<AGROW>
    DYSactive = [DYSactive subfreq(:,3,k)]; %#ok<AGROW>
    counterDYS = counterDYS+1;
end

Yrest(:,2) = mean(DYSrest,2);
% Erest(:,2)= std(DYSrest,0,2);
Erest(:,2) = std(DYSrest,0,2)/(sqrt(counterDYS));
Yactive(:,2) = mean(DYSactive,2);
% Eactive(:,2) = std(DYSactive,0,2);
Eactive(:,2) = std(DYSactive,0,2)/(sqrt(counterDYS));
%% grab ET data
% if ET data needs to be analyzed, change next line to "ET=true." Otherwise
% ET=false;
ET=true;
if ET
%     pathnameET = uigetdir('','Select directory with ET ecogPSD.mat files');
    etpath = uigetdir('', 'Select directory that contains ET _ecogPSD files to be analyzed');
    etpath = [etpath '\'];
    cd(etpath);
    
    ETdir = dir('*_ecogPSD.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script
    numET = length(ETdir);
    
    ETrest=[];
    ETactive=[];
    counterET=0;
    for i=1:numET
        filename = ETdir(i).name;
        load(filename);
        if isnan(order(k))
            continue; %skip files without contact pair over selected area
        end
        
     % Due to brain atrophy, sulci, vessels, etc, some subjects experieinced
% strong ECOG signals, while others had weaker signals. This creates a fair
% amount of variability in the data, given our small group sizes. To
% address this we will calculate the strength of the M1 signal relative to
% the parietal signal (contacts 1v2).
   for j = 1:size(rest.contact_pair,2)
    rest.contact_pair(1,j).rel_signal = rest.contact_pair(1,j).mean_PSD ./ ...
        mean(rest.contact_pair(1,1).mean_PSD);
    active.contact_pair(1,j).rel_signal = active.contact_pair(1,j).mean_PSD ...
        ./ mean(active.contact_pair(1,1).mean_PSD);
   end
    
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    [allfreq subfreq] = quantPSDshort_rel2parietal(rest,active,FREQ_QPSD,order,freq,filename);
    
    
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    %now create matrix of percent power in each freq band for each subject
    %to calculate "absolute" power in each band (rather than "% power") use
    %subfreq(:,1,k) and subfreq(:,3,k) for rest and active respectively
    ETrest = [ETrest subfreq(:,1,k)]; %#ok<AGROW>
    ETactive = [ETactive subfreq(:,3,k)]; %#ok<AGROW>
        counterET = counterET+1;
    end

    Yrest(:,3) = mean(ETrest,2);
%     Erest(:,3) = std(ETrest,0,2);
    Erest(:,3) = std(ETrest,0,2)/(sqrt(counterET));
    Yactive(:,3) = mean(ETactive,2);
%     Eactive(:,3) = std(ETactive,0,2);
    Eactive(:,3) = std(ETactive,0,2)/(sqrt(counterET));
end

%% grab OTHER data
%Added 4/6/09 for analysis of PD px with dystonic sx
% if OTHER data needs to be analyzed, change next line to "OTHER=true." Otherwise
% ET=false;
OTHER=false;
if OTHER
%     pathnameET = uigetdir('','Select directory with OTHER ecogPSD.mat files');
    OTHERpath = uigetdir('', 'Select directory that contains OTHER _ecogPSD files to be analyzed');
    OTHERpath = [OTHERpath '\'];
    cd(OTHERpath);
    
    OTHERdir = dir('*_ecogPSD.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script
    numOTHER = length(OTHERdir);
    
    OTHERrest=[];
    OTHERactive=[];
    counterOTHER=0;
    for i=1:numOTHER
        filename = OTHERdir(i).name;
        load(filename);
        if isnan(order(k))
            continue; %skip files without contact pair over selected area
        end
         % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
    [allfreq subfreq] = quantPSDshort(rest,active,FREQ_QPSD,order,freq,filename);
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    %now create matrix of percent power in each freq band for each subject
    OTHERrest = [OTHERrest subfreq(:,2,k)]; %#ok<AGROW>
        OTHERactive = [OTHERactive subfreq(:,4,k)]; %#ok<AGROW>
        counterOTHER = counterOTHER+1;
    end

    Yrest(:,3) = mean(OTHERrest,2);
%     Erest(:,3) = std(OTHERrest,0,2);
    Erest(:,3) = std(OTHERrest,0,2)/sqrt(counterOTHER);
    Yactive(:,3) = mean(OTHERactive,2);
%     Eactive(:,3) = std(OTHERactive,0,2);
    Eactive(:,3) = std(OTHERactive,0,2)/sqrt(counterOTHER);
end

%% plot data
hf=figure;
%temporary code for graphing only beta band
% Yrest= Yrest(2,:);
% Erest = Erest(2,:); 
% %------
subplot(2,1,1);
handlesr = barEB(Yrest,Erest);
% set(handlesr.ca,'XTick',[1;2;3;4;5],...
%     'XTickLabel',FREQBANDS,...
%     'ylim',YLIM,...
%     'YMinorGrid','on');
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'YMinorGrid','on');
title([BRAINAREAS{k} ' at REST']);
ylabel('Power');
hPD=handlesr.bars(1);
hDYS=handlesr.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
if ET
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['ET (n=' num2str(counterET) ')']);
    hET=handlesr.bars(3);
    set(hET,'FaceColor','m');
elseif OTHER
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['OTHER (n=' num2str(counterOTHER) ')']);
    hOTHER=handlesr.bars(3);
    set(hOTHER,'FaceColor','m');
else
%     legend(['PD (n=' num2str(counterPD) ')'],...
%         ['Dys (n=' num2str(counterDYS) ')']);
    legend(['PD'],...
        ['Dys']);
end

subplot(2,1,2);
handlesa=barEB(Yactive,Eactive);
% set(handlesa.ca,'XTick',[1;2;3;4;5],...
%     'XTickLabel',FREQBANDS,...
%     'ylim',YLIM,...
%     'YMinorGrid','on');
set(handlesa.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'YMinorGrid','on');
title([BRAINAREAS{k} ' during ' LIMBAREAS{l} ' MOVEMENT']);
ylabel('Power');
hPD=handlesa.bars(1);
hDYS=handlesa.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
if ET
    hET=handlesa.bars(3);
    set(hET,'FaceColor','m');
elseif OTHER
    hOTHER=handlesa.bars(3);
    set(hOTHER,'FaceColor','m');
end

function handles = barEB(barvals,errorvals)

% barEB plots bar graph of m-by-n matrix of bar values (by calling MATLAB
% function "bar") as well as error bars (by calling "errorbar")

% note: barvals and errorvals must be the same size!

% a list of handles are listed so that the user can change the properties
% of the plot
% handles.bars: handle to bar plot
% handles.gca: handle to axis

[numgroups numbars] = size(barvals); % # of groups, # bars in each group

handles.bars = bar(barvals);
hold on;
groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    handles.errors(i) = errorbar(x, barvals(:,i), errorvals(:,i), 'k', 'linestyle', 'none');
end

handles.ca = gca;
hold off;

%%
function [allfreq subfreq] = quantPSDshort_rel2parietal(rest,active,FREQ_QPSD,order,freq,filename)
allfreq = zeros(2,3,3);
for i=1:3
    if isnan(order(i))
        allfreq(:,:,i) = NaN;
        continue;
    end
%     rest_data = rest((order(i)),:); % from matrix-based code
%     active_data = active((order(i)),:);
    rest_data = rest.contact_pair(order(i)).rel_signal;
    active_data = active.contact_pair(order(i)).rel_signal;
% Need to limit rest and active data to those freq range of interest (typically 4Hz-100Hz)
%the _tmp variables are to limit the frequency spectrum for allfreq max
%values without affecting the later subfreq code, which also calls the
%rest, active, and freq variables

%to make more flexible, will use linear indexing to not require graphs to
%have 5 frequency bins
%     rest_data_tmp = rest_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
%     active_data_tmp = active_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
%     freq_tmp = freq(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));

%Need to change the _tmp such that the max freq is based on the freq bins
%and excludes, for example, 55-65Hz.
rest_data_tmp = rest_data(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));
active_data_tmp = active_data(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));
freq_tmp = freq(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));

    [c1 i1] = max(rest_data_tmp);
    allfreq(1,1,i)=freq_tmp(i1);
    allfreq(1,2,i)=c1;
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
for i = 1:size(FREQ_QPSD,1)
    subfreq = zeros(i,5,3);
end

for i=1:3
    if isnan(order(i))
    subfreq(:,:,i)=NaN;
    continue;
    end
    rest_data = rest.contact_pair(order(i)).rel_signal;
    active_data = active.contact_pair(order(i)).rel_signal;
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
% save('allfreq','subfreq');
% %% Save and write quantPSD data as excel spreadsheet
% [data4export] = exportdata(allfreq, subfreq, fn);
% cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.
return;
