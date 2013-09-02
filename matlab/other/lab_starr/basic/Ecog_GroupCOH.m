function Ecog_GroupCOH()
%The goal of this code is to import all transformed coherence files from
%each patient to create/plot average group coherence

%updated for the case that there is a valid m1 but not a valid s1 contact -
%now this situation will not stop the code--ALC

%ALC 12/18/09: updated to work with updated output from earlier analysis
%steps. Also updated to produce M1S1 coherence graphs from all 3 patient
%groups.

%% Set variables
FREQ_QCOH = [4 13;...     % delta alpha band
             13 22;...   % low beta band
             22 31;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   
% axis, xtick, and ytick parameters for plotting
AXIS = [0 50 0 1];
XTICK = 0:10:50;
YTICK = 0:.1:1;
YLIM = [0 .8];
AXIS_m1s1 = [0 100 0 1];
XTICK_m1s1 = 0:20:100;
%% Define frequency bands based on FREQ_QPSD 
%allows for user to create desired number of freq bins
FREQBANDS = [];
for i = 1:size(FREQ_QCOH,1)
    FREQBANDS = [FREQBANDS {num2str(FREQ_QCOH(i,:))}]; %#ok<AGROW>
end

%% Set directory path for finding relevant PD and Dys files and Import data
%This gets directory info for each file within the PD folder. For each file in directory, 
%need to import transcoh struct, select resting and active coh data under M1 contact 
%and place into separate (new) rest and active structures or matrices

LIMBAREAS = {'HAND' 'ELBOW' 'SHOULDER' 'JAW' 'FOOT' '"ARM"' '"NON-ARM"'};
l = menu('Select limb for analysis',LIMBAREAS);

% ------------PD-------------
% pdpath = uigetdir('', 'Select directory that contains PD _transcoh files to be analyzed');
% pdpath = [pdpath '\'];
% for 2 sec epoch, use line below
% pdpath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\PD';
% for full epoch length, use line below
pdpath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\PD\full epoch length';
pdpath = [pdpath '\' LIMBAREAS{l}];
cd(pdpath);

PDdir = dir('*_transcoh.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numPD = length(PDdir);

for i = 1:numPD
    PDcoh(i) = importdata(PDdir(i).name); %#ok<AGROW>
end

% ------------Dys--------------
% dyspath = uigetdir('', 'Select directory that contains Dys _transcoh files to be analyzed');
% dyspath = [dyspath '\'];
% for 2 sec epoch length, use line below
% dyspath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\DYS';
% for full length epoch length, use line below
dyspath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\DYS\full epoch length';
dyspath = [dyspath '\' LIMBAREAS{l}];
cd(dyspath);

DYSdir = dir('*_transcoh.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numDYS = length(DYSdir);

for i = 1:numDYS
    DYScoh(i) = importdata(DYSdir(i).name); %#ok<AGROW>
end
% -------------ET ------------- added 1/30/09
% etpath = uigetdir('', 'Select directory that contains ET _transcoh files to be analyzed');
% etpath = [etpath '\'];
% for 2 sec epoch, use line below:
% etpath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\ET';
% for full length epoch, use line below:
etpath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\ET\full epoch length';
etpath = [etpath '\' LIMBAREAS{l}];
cd(etpath);

ETdir = dir('*_transcoh.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numET = length(ETdir);

for i = 1:numET
    ETcoh(i) = importdata(ETdir(i).name); %#ok<AGROW>
end
%% Define S1 contact pairs

for i = 1:numPD
    PDcoh(i).S1 = PDcoh(i).M1-2; %#ok<AGROW>
end

for i = 1:numDYS
    DYScoh(i).S1 = DYScoh(i).M1-2; %#ok<AGROW>
end

for i = 1:numET
    ETcoh(i).S1 = ETcoh(i).M1-2; %#ok<AGROW>
end
%% Determine which PD subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numPD
    S1valuesPD(i) = PDcoh(i).S1; %#ok<AGROW> %defines S1 contact for each subject based on that subject's M1 contact
end
% A subject with contact 1 or 2 over M1 will not have an S1 contact, and
% must be excluded from S1 analysis
S1idx = find(S1valuesPD > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
S1values = nonzeros(S1valuesPD); %companion to S1idx: lists the S1 contacts for those patients who have valid S1 contacts
S1values = S1values'; 

%% Create matrix for aggregated PD S1-STN coherence data
PDrest_S1coh = [];
PDactive_S1coh = [];

for i = 1:length(S1idx)
    % Look for patients that have ecog-LFP coherence analysis
    if ~isfield(PDcoh(S1idx(i)).rest,'EcogLfp')
        continue;
    else
        PDrest_S1coh = [PDrest_S1coh PDcoh(S1idx(i)).rest.EcogLfp(:,int8(S1values(i)))]; %#ok<AGROW>
        PDactive_S1coh = [PDactive_S1coh PDcoh(S1idx(i)).active.EcogLfp(:,int8(S1values(i)))]; %#ok<AGROW>
    end
end
numPDs1 = size(PDrest_S1coh,2); %creates variable giving number of PD patients included in S1-STN analysis


%% Create matrix for aggregated PD M1-STN coherence data
PDrest_M1coh = []; 
PDactive_M1coh = [];

for i = 1: length(S1idx)
    % Look for patients that have ecog-LFP coherence analysis
    if ~isfield(PDcoh(S1idx(i)).rest,'EcogLfp')
        continue;
    else
        PDrest_M1coh(:,i) = PDcoh(i).rest.EcogLfp(:,int8(PDcoh(i).M1)); %#ok<AGROW> %obtains M1 coherence data
        PDactive_M1coh(:,i) = PDcoh(i).active.EcogLfp(:,int8(PDcoh(i).M1)); %#ok<AGROW>
    end
end
numPDm1 = size(PDrest_M1coh,2); %creates variable giving number of PD patients included in M1-STN analysis
%% Create matrix for aggregated PD M1-S1 coherence data
PDrest_M1S1 = [];
PDactive_M1S1 = [];

for i = 1:length(S1idx)
    % Look for patients that have M1-S1 coherence analysis
    if ~isfield(PDcoh(S1idx(i)).rest,'M1S1')
        continue;
    else
        PDrest_M1S1 = [PDrest_M1S1 PDcoh(S1idx(i)).rest.M1S1]; %#ok<AGROW>
        PDactive_M1S1 = [PDactive_M1S1 PDcoh(S1idx(i)).active.M1S1]; %#ok<AGROW>
    end
end
numPDm1s1 = size(PDrest_M1S1,2); %creates variable giving number of PD patients included in M1-S1 analysis

%% Determine which DYS subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numDYS
    S1valuesDYS(i) = DYScoh(i).S1; %#ok<AGROW> %defines S1 contact for each subject based on that subject's M1 contact
end
% A subject with contact 1 or 2 over M1 will not have an S1 contact, and
% must be excluded from S1 analysis
S1idx = find(S1valuesDYS > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
S1values = nonzeros(S1valuesDYS); %companion to S1idx: lists the S1 contacts for those patients who have valid S1 contacts
S1values = S1values'; 

%% Create matrix for aggregated DYS S1-STN coherence data
DYSrest_S1coh = [];
DYSactive_S1coh = [];

for i = 1:length(S1idx)
    % Look for patients that have ecog-LFP coherence analysis
    if ~isfield(DYScoh(S1idx(i)).rest,'EcogLfp')
        continue;
    else
        DYSrest_S1coh = [DYSrest_S1coh DYScoh(S1idx(i)).rest.EcogLfp(:,int8(S1values(i)))]; %#ok<AGROW>
        DYSactive_S1coh = [DYSactive_S1coh DYScoh(S1idx(i)).active.EcogLfp(:,int8(S1values(i)))]; %#ok<AGROW>
    end
end
numDYSs1 = size(DYSrest_S1coh,2);

%% Create matrix for aggregated DYS M1-STN coherence data
DYSrest_M1coh = []; 
DYSactive_M1coh = [];

for i = 1: numDYS
    % Look for patients that have ecog-LFP coherence analysis
    if ~isfield(DYScoh(S1idx(i)).rest,'EcogLfp')
        continue;
    else
        DYSrest_M1coh = [DYSrest_M1coh DYScoh(i).rest.EcogLfp(:,int8(DYScoh(i).M1))]; %#ok<AGROW>
        DYSactive_M1coh = [DYSactive_M1coh DYScoh(i).active.EcogLfp(:,int8(DYScoh(i).M1))]; %#ok<AGROW>
    end
end
numDYSm1 = size(DYSrest_M1coh,2);

%% Create matrix for aggregated DYS M1-S1 coherence data
DYSrest_M1S1 = [];
DYSactive_M1S1 = [];

for i = 1:length(S1idx)
    % Look for patients that have M1-S1 coherence analysis
    if ~isfield(DYScoh(S1idx(i)).rest,'M1S1')
        continue;
    else
        DYSrest_M1S1= [DYSrest_M1S1 DYScoh(S1idx(i)).rest.M1S1]; %#ok<AGROW>
        DYSactive_M1S1 = [DYSactive_M1S1 DYScoh(S1idx(i)).active.M1S1]; %#ok<AGROW>
    end
end
numDYSm1s1 = size(DYSrest_M1S1,2);

%% Determine which ET subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numET
    S1valuesET(i) = ETcoh(i).S1; %#ok<AGROW> %defines S1 contact for each subject based on that subject's M1 contact
end
% A subject with contact 1 or 2 over M1 will not have an S1 contact, and
% must be excluded from S1 analysis
S1idx = find(S1valuesET > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
S1values = nonzeros(S1valuesET); %companion to S1idx: lists the S1 contacts for those patients who have valid S1 contacts
S1values = S1values'; 

%% Create matrix for aggregated ET M1S1 coherence data
ETrest_M1S1 = [];
ETactive_M1S1 = [];

for i = 1:length(S1idx)
    % Look for patients that have M1-S1 coherence analysis
    if ~isfield(ETcoh(S1idx(i)).rest,'M1S1')
        continue;
    else
        ETrest_M1S1 = [ETrest_M1S1 ETcoh(S1idx(i)).rest.M1S1]; %#ok<AGROW>
        ETactive_M1S1 = [ETactive_M1S1 ETcoh(S1idx(i)).active.M1S1]; %#ok<AGROW>
    end
end
numETm1s1 = size(ETrest_M1S1,2);


%% Average across patients
% ------------PD-----------------
Mean_PDrest_S1coh = mean(PDrest_S1coh,2);
Mean_PDactive_S1coh = mean(PDactive_S1coh,2);
Mean_PDrest_M1coh = mean(PDrest_M1coh,2);
Mean_PDactive_M1coh = mean(PDactive_M1coh,2);
Mean_PDrest_M1S1 = mean(PDrest_M1S1,2);
Mean_PDactive_M1S1 = mean(PDactive_M1S1,2);
% ----------Dystonia----------------
Mean_DYSrest_S1coh = mean(DYSrest_S1coh,2);
Mean_DYSactive_S1coh = mean(DYSactive_S1coh,2);
Mean_DYSrest_M1coh = mean(DYSrest_M1coh,2);
Mean_DYSactive_M1coh = mean(DYSactive_M1coh,2);
Mean_DYSrest_M1S1 = mean(DYSrest_M1S1,2);
Mean_DYSactive_M1S1 = mean(DYSactive_M1S1,2);
% ------------ET-----------------
Mean_ETrest_M1S1 = mean(ETrest_M1S1,2);
Mean_ETactive_M1S1 = mean(ETactive_M1S1,2);
% %To plot frequency, it must be in array format, not structure format
F = PDcoh.freq; 

%% Quantify transformed coherence levels in each frequency band using quantCOH
% quantCOH takes aggregated coherence data and calculates population mean
% and standard error of all frequency bands
%--------M1-STN coherence------------
Yrest_M1 = zeros(size(FREQ_QCOH,1),2);
Erest_M1 = zeros(size(FREQ_QCOH,1),2);

[Yrest_M1(:,1) Erest_M1(:,1)] = quantCOH(PDrest_M1coh,FREQ_QCOH,F);
[Yrest_M1(:,2) Erest_M1(:,2)] = quantCOH(DYSrest_M1coh,FREQ_QCOH,F);

[Yactive_M1(:,1) Eactive_M1(:,1)] = quantCOH(PDactive_M1coh,FREQ_QCOH,F);
[Yactive_M1(:,2) Eactive_M1(:,2)] = quantCOH(DYSactive_M1coh,FREQ_QCOH,F);

%--------S1-STN coherence------------
Yrest_S1 = zeros(size(FREQ_QCOH,1),2);
Erest_S1 = zeros(size(FREQ_QCOH,1),2);

[Yrest_S1(:,1) Erest_S1(:,1)] = quantCOH(PDrest_S1coh,FREQ_QCOH,F);
[Yrest_S1(:,2) Erest_S1(:,2)] = quantCOH(DYSrest_S1coh,FREQ_QCOH,F);

[Yactive_S1(:,1) Eactive_S1(:,1)] = quantCOH(PDactive_S1coh,FREQ_QCOH,F);
[Yactive_S1(:,2) Eactive_S1(:,2)] = quantCOH(DYSactive_S1coh,FREQ_QCOH,F);

%--------M1-S1 coherence------------
Yrest_M1S1 = zeros(size(FREQ_QCOH,1),3);
Erest_M1S1 = zeros(size(FREQ_QCOH,1),3);

PDrestM1S1 = quantFREQ(PDrest_M1S1,FREQ_QCOH,F);
DYSrestM1S1 = quantFREQ(DYSrest_M1S1,FREQ_QCOH,F);
ETrestM1S1 = quantFREQ(ETrest_M1S1,FREQ_QCOH,F);
[Yrest_M1S1(:,1) Erest_M1S1(:,1)] = quantCOH(PDrest_M1S1,FREQ_QCOH,F);
[Yrest_M1S1(:,2) Erest_M1S1(:,2)] = quantCOH(DYSrest_M1S1,FREQ_QCOH,F);
[Yrest_M1S1(:,3) Erest_M1S1(:,3)] = quantCOH(ETrest_M1S1,FREQ_QCOH,F);


PDactiveM1S1 = quantFREQ(PDactive_M1S1,FREQ_QCOH,F);
DYSactiveM1S1 = quantFREQ(DYSactive_M1S1,FREQ_QCOH,F);
ETactiveM1S1 = quantFREQ(ETactive_M1S1,FREQ_QCOH,F);
[Yactive_M1S1(:,1) Eactive_M1S1(:,1)] = quantCOH(PDactive_M1S1,FREQ_QCOH,F);
[Yactive_M1S1(:,2) Eactive_M1S1(:,2)] = quantCOH(DYSactive_M1S1,FREQ_QCOH,F);
[Yactive_M1S1(:,3) Eactive_M1S1(:,3)] = quantCOH(ETactive_M1S1,FREQ_QCOH,F);


%% Plot M1 only
hf1 = figure;
subplot(2,1,1);
plot(F, Mean_PDrest_M1coh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSrest_M1coh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
axis(AXIS) 
set (gca,'xtick',XTICK)
set (gca,'ytick',YTICK)
ylabel('Transformed Coherence')
xlabel('F')
legend(['PD (n=',num2str(numPDm1),')'], ['DYS (n=' num2str(numDYSm1) ')']);
title ('Transformed coherence, M1-STN during REST')
    
subplot(2,1,2);
plot(F, Mean_PDactive_M1coh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSactive_M1coh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
axis(AXIS) 
set(gca,'xtick',XTICK)
set(gca,'ytick',YTICK)
ylabel('Transformed Coherence')
xlabel('F')
legend(['PD (n=',num2str(numPDm1),')'], ['DYS (n=' num2str(numDYSm1) ')']);
title (['Transformed coherence, M1-STN during ' LIMBAREAS{l} ' MOVEMENT']);
%% Second Plot for M1 only
hf2 = figure;

subplot(2,1,1);
plot(F, Mean_PDrest_M1coh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_PDactive_M1coh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
axis(AXIS);
set (gca,'xtick',XTICK);
set (gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['REST (n=',num2str(numPDm1),')'], ['ACTIVE ' LIMBAREAS{l}]);
title ('Transformed coherence, M1-STN in PD');
    
subplot(2,1,2);
plot(F, Mean_DYSrest_M1coh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSactive_M1coh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
axis(AXIS);
set(gca,'xtick',XTICK);
set(gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['REST (n=',num2str(numDYSm1),')'], ['ACTIVE ' LIMBAREAS{l}]);
title ('Transformed coherence, M1-STN in Dystonia');
    
%% Plot S1 only
hf3 = figure;
subplot(2,1,1);
plot(F, Mean_PDrest_S1coh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSrest_S1coh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
axis(AXIS);
set(gca,'xtick',XTICK);
set(gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['PD (n=',num2str(numPDs1),')'], ['DYS (n=' num2str(numDYSs1) ')']);
title ('Transformed coherence, S1-STN during REST');

subplot(2,1,2);
plot(F, Mean_PDactive_S1coh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSactive_S1coh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
axis(AXIS);
set (gca,'xtick',XTICK);
set (gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['PD (n=',num2str(numPDs1),')'], ['DYS (n=' num2str(numDYSs1) ')']);
title (['Transformed coherence, S1-STN during ' LIMBAREAS{l} ' MOVEMENT']);
%% Second Plot for S1 only
hf4 = figure;

subplot(2,1,1);
plot(F, Mean_PDrest_S1coh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_PDactive_S1coh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
axis(AXIS);
set (gca,'xtick',XTICK);
set (gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['REST (n=',num2str(numPDs1),')'], ['ACTIVE ' LIMBAREAS{l}]);
title ('Transformed coherence, S1-STN in PD');

subplot(2,1,2);
plot(F, Mean_DYSrest_S1coh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSactive_S1coh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
axis(AXIS);
set (gca,'xtick',XTICK);
set (gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['REST (n=',num2str(numDYSs1),')'],['ACTIVE ' LIMBAREAS{l}]);
title ('Transformed coherence, S1-STN in Dystonia');
    
%% Plot M1S1 coherence 
% Plot M1S1 #1
hf5 = figure; 
subplot(2,1,1);
plot(F, Mean_PDrest_M1S1, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSrest_M1S1, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
hold on;
plot(F, Mean_ETrest_M1S1, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
axis(AXIS_m1s1)
set (gca,'xtick',XTICK_m1s1)
set (gca,'ytick',YTICK)
ylabel('Transformed Coherence')
xlabel('F')
legend(['PD (n=',num2str(numPDm1s1),')'], ['DYS (n=' num2str(numDYSm1s1) ')'],...
    ['ET (n=',num2str(numETm1s1),')']);
title ('Transformed coherence, M1-S1 during REST')
    
subplot(2,1,2);
plot(F, Mean_PDactive_M1S1, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSactive_M1S1, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
hold on;
plot(F, Mean_ETactive_M1S1, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
axis(AXIS_m1s1);
set (gca,'xtick',XTICK_m1s1);
set (gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['PD (n=',num2str(numPDm1s1),')'], ['DYS (n=' num2str(numDYSm1s1) ')'],...
    ['ET (n=',num2str(numETm1s1),')']);
title (['Transformed coherence, M1-S1 during ' LIMBAREAS{l} ' MOVEMENT']);

%% Plot #2 for M1S1
hf6 = figure;

subplot(3,1,1);
plot(F, Mean_PDrest_M1S1, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_PDactive_M1S1, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
axis(AXIS_m1s1);
set (gca,'xtick',XTICK_m1s1);
set (gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['REST (n=',num2str(numPDm1s1),')'], ['ACTIVE ' LIMBAREAS{l}]);
title ('Transformed coherence, M1-S1 in PD');
    
subplot(3,1,2);
plot(F, Mean_DYSrest_M1S1, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSactive_M1S1, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
axis(AXIS_m1s1);
set (gca,'xtick',XTICK_m1s1);
set (gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['REST (n=',num2str(numDYSm1s1),')'], ['ACTIVE ' LIMBAREAS{l}]);
title ('Transformed coherence, M1-S1 in Dystonia');
    
subplot(3,1,3);
plot(F, Mean_ETrest_M1S1, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_ETactive_M1S1, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','r','LineWidth',1.5);
axis(AXIS_m1s1);
set (gca,'xtick',XTICK_m1s1);
set (gca,'ytick',YTICK);
ylabel('Transformed Coherence');
xlabel('F');
legend(['REST (n=',num2str(numETm1s1),')'], ['ACTIVE ' LIMBAREAS{l}]);
title ('Transformed coherence, M1-S1 in ET');

%% plot quantCOH results

%--------M1-STN coherence------------
hf7=figure;
subplot(2,1,1);
handlesr = barEB(Yrest_M1,Erest_M1);
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YTick',YTICK,...
    'YMinorGrid','on');
title('Transformed coherence, M1-STN at REST');
ylabel('Transformed coherence');
hPD=handlesr.bars(1);
hDYS=handlesr.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
legend(['PD (n=' num2str(numPDm1) ')'],...
        ['Dys (n=' num2str(numDYSm1) ')']);

subplot(2,1,2);
handlesa=barEB(Yactive_M1,Eactive_M1);
set(handlesa.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YTick',YTICK,...
    'YMinorGrid','on');
title(['Transformed coherence, M1-STN during ' LIMBAREAS{l} ' MOVEMENT']);
ylabel('Transformed coherence');
hPD=handlesa.bars(1);
hDYS=handlesa.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');

%--------S1-STN coherence------------
hf8=figure;
subplot(2,1,1);
handlesr = barEB(Yrest_S1,Erest_S1);
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YTick',YTICK,...
    'YMinorGrid','on');
title('Transformed coherence, S1-STN at REST');
ylabel('Transformed coherence');
hPD=handlesr.bars(1);
hDYS=handlesr.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
legend(['PD (n=' num2str(numPDs1) ')'],...
        ['Dys (n=' num2str(numDYSs1) ')']);

subplot(2,1,2);
handlesa=barEB(Yactive_S1,Eactive_S1);
set(handlesa.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YTick',YTICK,...
    'YMinorGrid','on');
title(['Transformed coherence, S1-STN during ' LIMBAREAS{l} ' MOVEMENT']);
ylabel('Transformed coherence');
hPD=handlesa.bars(1);
hDYS=handlesa.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');

%--------M1-S1 coherence------------
hf9=figure;
subplot(2,1,1);
handlesr = barEB(Yrest_M1S1,Erest_M1S1);
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YTick',YTICK,...
    'YMinorGrid','on');
title('Transformed coherence, M1-S1 at REST');
ylabel('Transformed coherence');
hPD=handlesr.bars(1);
hDYS=handlesr.bars(2);
hET=handlesr.bars(3);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
set(hET,'FaceColor','m');
legend(['PD (n=' num2str(numPDm1s1) ')'],...
        ['Dys (n=' num2str(numDYSm1s1) ')'],...
        ['ET (n=' num2str(numETm1s1) ')']);

subplot(2,1,2);
handlesa=barEB(Yactive_M1S1,Eactive_M1S1);
set(handlesa.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YTick',YTICK,...
    'YMinorGrid','on');
title(['Transformed coherence, M1-S1 during ' LIMBAREAS{l} ' MOVEMENT']);
ylabel('Transformed coherence');
hPD=handlesa.bars(1);
hDYS=handlesa.bars(2);
hET=handlesa.bars(3);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
set(hET,'FaceColor','m');

%% save variables
cd('C:\Users\Starr\Documents\ECOG data\Trans_coh_data\ANOVA stats');
save(['Ecog_GroupCOH_' LIMBAREAS{l}],...
    'Yrest_M1','Erest_M1','Yactive_M1','Eactive_M1',...
    'Yrest_S1','Erest_S1','Yactive_S1','Eactive_S1',...
    'PDrestM1S1','PDactiveM1S1','DYSrestM1S1','DYSactiveM1S1',...
    'ETrestM1S1','ETactiveM1S1',...
    'Yrest_M1S1','Erest_M1S1','Yactive_M1S1','Eactive_M1S1');
%%
% %% Plot Premotor, M1, and M1S1
% hf = figure;
% 
% for i = 1:3
% subplot (2,3,i);
% plot(F, Mean_PD_rest_tcoh(:,:,i), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_DYS_rest_tcoh(:,:,i), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 50 0 1]) %axis of alpha and beta freq
% %     set (gca,'xtick',[0:20:150])
% %     set (gca,'xtick',[0:20:150])
% %     set (gca,'ytick',[0:0.2:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
%     if i == 1
%         title('Transformed coherence over Premotor Cortex')
%     elseif i == 2
%         title ('Transformed coherece over M1')
%     elseif i == 3
%         title ('Transformed coherece over S1-M1')
%     end
%     
%     
% subplot(2,3,i+3);
% plot(F, Mean_PD_active_tcoh(:,:,i), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_DYS_active_tcoh(:,:,i), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 50 0 1]) %axis of alpha and beta freq
% %     set (gca,'xtick',[0:20:150])
% %     set (gca,'xtick',[0:20:150])
% %     set (gca,'ytick',[0:0.2:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
% %     title(['Transformed coherence during Movement: PD vs DYS'])
% end
% 
% % % Create textbox
% % annotation(figure1,'textbox','String',{'Rest','Condition'},...
% %     'FitHeightToText','off',...
% %     'Position',[0.004171 0.7504 0.07821 0.06317]);
% % 
% % % Create textbox
% % annotation(figure1,'textbox','String',{'Active','Movement'},...
% %     'FitHeightToText','off',...
% %     'Position',[0.002086 0.2496 0.07404 0.06317]);

% %% Second plot Premotor, M1, and M1S1
% hf2 = figure;
% 
% for i = 1:3
% subplot (2,3,i);
% plot(F, Mean_PD_rest_tcoh(:,:,i), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_PD_active_tcoh(:,:,i), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 50 0 1]) %axis of alpha and beta freq
% %     set (gca,'xtick',[0:20:150])
% %     set (gca,'xtick',[0:20:150])
% %     set (gca,'ytick',[0:0.2:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['PD rest (n=',num2str(numPD),')'], ['PD active']);
%     if i == 1
%         title('Transformed coherence over Premotor Cortex')
%     elseif i == 2
%         title ('Transformed coherece over M1')
%     elseif i == 3
%         title ('Transformed coherece over S1-M1')
%     end
%     
%     
% subplot(2,3,i+3);
% plot(F, Mean_DYS_rest_tcoh(:,:,i), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_DYS_active_tcoh(:,:,i), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 50 0 1]) %axis of alpha and beta freq
% %     set (gca,'xtick',[0:20:150])
% %     set (gca,'xtick',[0:20:150])
% %     set (gca,'ytick',[0:0.2:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['Dys rest (n=',num2str(numDYS),')'], ['Dys active']);
% %     title(['Transformed coherence during Movement: PD vs DYS'])
% end
% 
% % % Create textbox
% % annotation(figure1,'textbox','String',{'PD'},...
% %     'FitHeightToText','off',...
% %     'Position',[0.004171 0.7504 0.07821 0.06317]);
% % 
% % % Create textbox
% % annotation(figure1,'textbox','String',{'Dys'},...
% %     'FitHeightToText','off',...
% %     'Position',[0.002086 0.2496 0.07404 0.06317]);

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

% return;
