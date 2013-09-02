%The goal of this code is to import all transformed coherence files from
%each patient to create/plot average group coherence

%% Set directory path for finding relevant PD and Dys files and Import data
%This gets directory info for each file within the PD folder. For each file in directory, 
%need to import transcoh struct, select resting and active coh data under M1 contact 
%and place into separate (new) rest and active structures or matrices

% for PD
pdpath = uigetdir('', 'Select directory that contains PD _transcoh files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numPD = length(PDdir);

for i = 1:numPD
    PDcoh(i) = importdata(PDdir(i).name);
end

% for Dys
dyspath = uigetdir('', 'Select directory that contains Dys _transcoh files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numDYS = length(DYSdir);

for i = 1:numDYS
    DYScoh(i) = importdata(DYSdir(i).name);
end

%% Aggregate all M1 coherence data from Rest state and Active state

% for PD
PD_rest_M1 = zeros(length(PDcoh(1).rest),numPD);
PD_active_M1 = zeros(length(PDcoh(1).active),numPD);
for i = 1: numPD
    PD_rest_M1(:,i) = PDcoh(i).rest(:,M1_contact);
end

for i = 1:numPD
    PD_active_M1(:,i) = PDcoh(i).active(:,M1_contact);
end

% for Dys
DYS_rest_M1 = zeros(length(DYScoh(1).rest),numDYS);
DYS_active_M1 = zeros(length(DYScoh(1).active),numDYS);

for i = 1:numDYS
    DYS_rest_M1(:,i) = DYScoh(i).rest(:,M1_contact);
end

for i = 1:numDYS
    DYS_active_M1(:,i) = DYScoh(i).active(:,M1_contact);
end

%% Average across patients
Mean_PD_rest_tcoh = mean(PD_rest_M1,2);
Mean_PD_active_tcoh = mean(PD_active_M1,2);
Mean_DYS_rest_tcoh = mean(DYS_rest_M1,2);
Mean_DYS_active_tcoh = mean(DYS_active_M1,2);

% %To plot frequency, it must be in array format, not structure format
F = PDcoh.freq; 
%% Plot
hf = figure;

subplot (2,1,1);
plot(F, Mean_PD_rest_tcoh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_rest_tcoh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
    axis([0 150 0 1]) % wide view axis
%     axis([0 50 0 1]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
%     set (gca,'xtick',[0:20:150])
%     set (gca,'ytick',[0:0.2:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
    title(['Transformed coherence at Rest: PD vs DYS'])
    
subplot(2,1,2);
plot(F, Mean_PD_active_tcoh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_active_tcoh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
    axis([0 150 0 1]) % wide view axis
%     axis([0 50 0 1]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
%     set (gca,'xtick',[0:20:150])
%     set (gca,'ytick',[0:0.2:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
    title(['Transformed coherence during Movement: PD vs DYS'])