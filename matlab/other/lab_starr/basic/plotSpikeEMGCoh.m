function plotSpikeEMGCoh()
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes in names of excel file/sheet that contains the data 
% output from Spike-EMG coherence analysis (RunSpikeEMGCoh.m), then outputs
% a graphs of frequency distribution and text containing relevant data.
% 
% CREATED BY: Sho Shimamoto
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define important parameters
EDGES = [1e-3 1 2 3 5 10 15 20 30 60 100 200]; 
BIN_RANGE = {'0-1','1-2','2-3','3-5','5-10','10-15','15-20','20-30','30-60','60-100','100-200 ',''};
FILE = {'SpikeEMGCohGPi.xls' 'SpikeEMGCohGPe.xls' 'SpikeEMGCohSTN.xls' 'SpikeEMGCohTH.xls'};
CH = {'Ch 1' 'Ch 2' 'Ch 3' 'Ch 4' 'Ch 5'};

% change directory to one containing SpikeEMGCoh excel files
cd('C:\Documents and Settings\HP_Administrator\My Documents\Lab Documents\Data\SpreadSHEETs\SpikeEMGCoh');
% note: another option is to allow user to manually select directory by
% commenting out the command cd above and uncommenting the next line:
% cd(uigetdir);

% Create the str 'filename'
filename = FILE{menu('Select file', FILE)};

% Create SHEETname options for the menu function
if strcmp(filename,'SpikeEMGCohGPi.xls')
    SHEET = {'Dyst GPi spont' 'Dyst GPi active' 'Dyst GPi task'...
        'PD GPi spont' 'PD GPi active'...
        'Chorea GPi spont' 'Chorea GPi active'...
        'Midbrain tremor GPi spont' 'Midbrain tremor GPi active'...
        'nNHP GPi spont' 'nNHP GPi task' 'dysNHP GPi spont' 'dysNHP GPi task'};
elseif strcmp(filename,'SpikeEMGCohGPe.xls')
    SHEET = {'Dyst GPe spont'...
        'PD GPe spont' 'Chorea GPe spont'...
        'Midbrain tremor GPe spont' 'nNHP GPe spont' 'nNHP GPe task'};
elseif strcmp(filename,'SpikeEMGCohSTN.xls')
    SHEET = {'Dyst STN spont' 'Dyst STN active'...
        'PD STN spont' 'PD STN active' 'PD STN task'};
elseif strcmp(filename, 'SpikeEMGCohTH.xls')
    SHEET = {'Midbrain tremor TH spont' 'Midbrain tremor TH active'};
end

% Create the str 'SHEETname'
sheetname = SHEET{menu('Select SHEET', SHEET)};


% Ask which channel should have its oscillation eliminated
% note: Ch 5 (accel ch) is automatically eliminated
ch_elim = menu('Select channel to eliminate from Coh data', CH);

% Import Spike-EMG coh data
[data, type] = xlsread(filename, sheetname);


%% spike auto-correlation
%frequency distribution
[numunits,numcols] = size(data); % calculate # units being analyzed and # columns
numsigunits = sum(data(:,2)>0);  % calculate number of units with at least one significant AC
percentsigunits = numsigunits/numunits * 100; % calculate percent of units with significant AC
tmp = [data(:,2) data(:,5)];
AC.freq = tmp(tmp>0);
n = 100*histc(AC.freq, EDGES)/numunits;

% crate histogram
figure1 = figure;
axes('Parent', figure1,...
    'XTickLabel',BIN_RANGE,...
    'XTick', 1:12);
box('on')
hold('all')
bar(n);
xlim([0 12]);
ylim('auto');
xlabel('bins (Hz)'); ylabel ('Prop of units studied (%)','fontsize',10);
title([sheetname, '  Auto-Correlation','  # units=', int2str(numunits),  ', prop sig oscil=', int2str(percentsigunits),'%']);

% stats
AC.mean = mean(AC.freq);
AC.median = median(AC.freq);
AC.std = std(AC.freq);
%% Spike-EMG cross-correlation
% frequency distribution
% count total units studied w/ EMG
numEMGunits = sum(data(:,9)>0);
% note:eliminate accel Ch 5 and ch_elim
counter = 10;
tmp = [];
XC.freq = [];
while counter < numcols
    tmp = [tmp data(:,counter)];
    XC.freq = [XC.freq data(:,counter+1)];
    counter = counter + 4;
end

XC.freq = XC.freq(tmp~=5 & tmp~=ch_elim); 
XC.freq = XC.freq(XC.freq>0);

% crate histogram
figure2 = figure;
n = histc(XC.freq, EDGES);
axes('Parent', figure2,...
    'XTickLabel',BIN_RANGE,...
    'XTick', 1:12);
box('on')
hold('all')
bar(n);
xlim([0 12]);
ylim('auto');
xlabel('bins (Hz)'); ylabel ('# units','fontsize',10);
title([sheetname, '  Cross-Correlation','  # units=', int2str(numEMGunits)]);

% proportion of units w/ at least one significant coh
prop_1 = data(:,11)/numEMGunits;

% also, give proportion of units w/ coh in 4, 3, 2, or just 1 channel

% note: store data in XC so info can be accessed later





