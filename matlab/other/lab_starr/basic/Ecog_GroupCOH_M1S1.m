function Ecog_GroupCOH_M1S1()
%The goal of this code is to import all transformed coherence files from
%each patient to create/plot average group coherence

%updated for the case that there is a valid m1 but not a valid s1 contact -
%now this situation will not stop the code--ALC

% Created by SAS 5/5/2010 to plot M1-S1 Coherence ONLY (ie. no M1-LFP or
% S1-LFP coherence).

%% Set directory path for finding relevant PD and Dys files and Import data
%This gets directory info for each file within the PD folder. For each file in directory, 
%need to import transcoh struct, select resting and active coh data under M1 contact 
%and place into separate (new) rest and active structures or matrices

LIMBAREAS = {'HAND' 'ELBOW' 'SHOULDER' 'JAW' 'FOOT' '"ARM"' '"NON-ARM"'};
l = menu('Select limb for analysis',LIMBAREAS);

% for PD
% pdpath = uigetdir('', 'Select directory that contains PD _transcoh files to be analyzed');
% pdpath = [pdpath '\'];
pdpath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\PD';
pdpath = [pdpath '\' LIMBAREAS{l}];
cd(pdpath);

PDdir = dir('*_transcoh.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numPD = length(PDdir);

for i = 1:numPD
    PDcoh(i) = importdata(PDdir(i).name);
end

% for Dys
% dyspath = uigetdir('', 'Select directory that contains Dys _transcoh files to be analyzed');
% dyspath = [dyspath '\'];
dyspath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\DYS';
dyspath = [dyspath '\' LIMBAREAS{l}];
cd(dyspath);

DYSdir = dir('*_transcoh.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numDYS = length(DYSdir);

for i = 1:numDYS
    DYScoh(i) = importdata(DYSdir(i).name);
end
% for ET - added 1/30/09
% etpath = uigetdir('', 'Select directory that contains ET _transcoh files to be analyzed');
% etpath = [etpath '\'];
etpath = 'C:\Users\Starr\Documents\ECOG data\Trans_coh_data\ET';
etpath = [etpath '\' LIMBAREAS{l}];
cd(etpath);

ETdir = dir('*.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numET = length(ETdir);

for i = 1:numET
    ETcoh(i) = importdata(ETdir(i).name);
end
%% Define S1 contact pairs

for i = 1:numPD
%     PDcoh(i).S1 = int8(PDcoh(i).rest_M1) - 2;
    PDcoh(i).S1 = PDcoh(i).M1-2;
end

for i = 1:numDYS
%     DYScoh(i).S1 = int8(DYScoh(i).rest_M1) - 2;
    DYScoh(i).S1 = DYScoh(i).M1-2;
end

for i = 1:numET
%     ETcoh(i).S1 = int8(ETcoh(i).rest_M1) - 2;
    ETcoh(i).S1 = ETcoh(i).M1-2;
end
%% Determine which PD subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numPD
    S1valuesPD(i) = PDcoh(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
end
% A subject with contact 1 or 2 over M1 will not have an S1 contact, and
% must be excluded from S1 analysis
S1idx = find(S1valuesPD > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
S1values = nonzeros(S1valuesPD); %companion to S1idx: lists the S1 contacts for those patients who have valid S1 contacts
S1values = S1values'; 
numPDs1 = length(S1idx);

%% Create matrix for aggregated PD M1S1 coherence data
PDrest_M1S1 = zeros(length(PDcoh(1).rest.M1S1),numPDs1); 
PDactive_M1S1 = zeros(length(PDcoh(1).active.M1S1),numPDs1);

for i = 1:length(S1idx)
    PDrest_M1S1(:,i) = PDcoh(S1idx(i)).rest.M1S1;
    PDactive_M1S1(:,i) = PDcoh(S1idx(i)).active.M1S1;
end

%% Determine which DYS subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numDYS
    S1valuesDYS(i) = DYScoh(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
end
% A subject with contact 1 or 2 over M1 will not have an S1 contact, and
% must be excluded from S1 analysis
S1idx = find(S1valuesDYS > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
S1values = nonzeros(S1valuesDYS); %companion to S1idx: lists the S1 contacts for those patients who have valid S1 contacts
S1values = S1values'; 
numDYSs1=length(S1idx);


%% Create matrix for aggregated DYS M1S1 coherence data
DYSrest_M1S1 = zeros(length(DYScoh(1).rest.M1S1), numDYSs1);
DYSactive_M1S1 = zeros(length(DYScoh(1).active.M1S1),numDYSs1);

for i = 1:length(S1idx)
    DYSrest_M1S1(:,i) = DYScoh(S1idx(i)).rest.M1S1;
    DYSactive_M1S1(:,i) = DYScoh(S1idx(i)).active.M1S1;
end
%% Determine which ET subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numET
    S1valuesET(i) = ETcoh(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
end
% A subject with contact 1 or 2 over M1 will not have an S1 contact, and
% must be excluded from S1 analysis
S1idx = find(S1valuesET > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
S1values = nonzeros(S1valuesET); %companion to S1idx: lists the S1 contacts for those patients who have valid S1 contacts
S1values = S1values'; 
numETs1=length(S1idx);

%% Create matrix for aggregated ET M1S1 coherence data
ETrest_M1S1 = zeros(length(ETcoh(1).rest.M1S1), numETs1);
ETactive_M1S1 = zeros(length(ETcoh(1).active.M1S1), numETs1);

for i = 1:length(S1idx)
    ETrest_M1S1(:,i) = ETcoh(S1idx(i)).rest.M1S1;
    ETactive_M1S1(:,i) = ETcoh(S1idx(i)).active.M1S1;
end
numETs1 = length(ETrest_M1S1(1,:));

% %% Create matrix for aggregated ET S1 coherence data
% ETrest_S1coh = zeros(length(ETcoh(1).rest.EcogLfp), length(S1idx));
% ETactive_S1coh = zeros(length(ETcoh(1).active.EcogLfp), length(S1idx));
% 
% for i = 1:length(S1idx)
%     ETrest_S1coh(:,i) = ETcoh(S1idx(i)).rest.EcogLfp(:,int8(S1values(i)));
%     ETactive_S1coh(:,i) = ETcoh(S1idx(i)).active.EcogLfp(:,int8(S1values(i)));
% end
% numETs1 = length(ETrest_S1coh(1,:));
% %% Create matrix for aggregated ET M1 coherence data
% ETrest_M1coh = zeros(length(ETcoh(1).rest.EcogLfp),numET); 
% ETactive_M1coh = zeros(length(ETcoh(1).active.EcogLfp),numET);
% 
% for i = 1:numET
% ETrest_M1coh(:,i) = ETcoh(i).rest.EcogLfp(:,int8(ETcoh(i).M1)); %obtains M1 coherence data
% ETactive_M1coh(:,i) = ETcoh(i).active.EcogLfp(:,int8(ETcoh(i).M1));
% end

%% Average across patients
Mean_PDrest_M1S1 = mean(PDrest_M1S1,2);
Mean_PDactive_M1S1 = mean(PDactive_M1S1,2);
Mean_DYSrest_M1S1 = mean(DYSrest_M1S1,2);
Mean_DYSactive_M1S1 = mean(DYSactive_M1S1,2);
Mean_ETrest_M1S1 = mean(ETrest_M1S1,2);
Mean_ETactive_M1S1 = mean(ETactive_M1S1,2);
% %To plot frequency, it must be in array format, not structure format
F = PDcoh.freq; 


%% Plot M1S1 coherence 
% Plot M1S1 #1
hf5 = figure; 
subplot(2,1,1);
plot(F, Mean_PDrest_M1S1, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSrest_M1S1, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
hold on;
plot(F, Mean_ETrest_M1S1, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
    axis([0 100 0 1]) 
    set (gca,'xtick',[0:20:100])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')'],...
        ['ET (n=',num2str(numET),')']);
%     legend(['PD'], ['DYS']);
    title ('Transformed M1S1 coherence at REST');
    
subplot(2,1,2);
plot(F, Mean_PDactive_M1S1, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSactive_M1S1, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
hold on;
plot(F, Mean_ETactive_M1S1, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
    axis([0 100 0 1]) 
    set (gca,'xtick',[0:20:100])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')'],...
        ['ET (n=',num2str(numET),')']);
%     legend(['PD'], ['DYS']);
    title (['Transformed M1-S1 coherence during ' LIMBAREAS{l} ' MOVEMENT']);

% Plot #2 for M1S1
hf6 = figure;

subplot(3,1,1);
plot(F, Mean_PDrest_M1S1, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_PDactive_M1S1, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
    axis([0 100 0 1]) 
    set (gca,'xtick',[0:20:100])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD rest (n=',num2str(numPDs1),')'], ['PD active']);
    title (['Transformed M1-S1 coherence in PD, ' LIMBAREAS{l} ' MOVEMENT']);
    
subplot(3,1,2);
plot(F, Mean_DYSrest_M1S1, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYSactive_M1S1, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
    axis([0 100 0 1]) ;
    set (gca,'xtick',[0:20:100]);
    set (gca,'ytick',[0:.1:1]);
    ylabel('Transformed Coherence');
    xlabel('F');
    legend(['Dys rest (n=',num2str(numDYSs1),')'], ['Dys active']);
    title ('Transformed M1S1 coherence in Dystonia');
    
subplot(3,1,3);
plot(F, Mean_ETrest_M1S1, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_ETactive_M1S1, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','r','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
    axis([0 100 0 1]);
    set (gca,'xtick',[0:20:100]);
    set (gca,'ytick',[0:.1:1]);
    ylabel('Transformed Coherence');
    xlabel('F');
    legend(['ET rest (n=',num2str(numETs1),')'], ['ET active']);
    title ('Transformed M1S1 coherence in ET');
