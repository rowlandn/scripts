%The goal of this code is to import all transformed coherence files from
%each patient to create/plot average group coherence
%Originally called 'groupM1S1coh.m'

%UPDATE: ALC 6/11/09: name changed to M1S1_GroupCOH.mat and code updated to
%match the updates made in input files. Specifically this code now looks
%for *_transcoh.mat files rather than *_m1s1coh.mat files. The name of the
%variable where M1S1 data is stored by upstream data analysis has changed,
%so this group analysis code has been updated be able to handle the new
%variable name. 

%INPUT: *_transcoh.mat
%OUTPUT: graphs

%% Menu of limb areas that can be used in analysis
LIMBAREAS = {'HAND' 'ELBOW' 'SHOULDER' 'JAW' 'FOOT' '"ARM"' '"NON-ARM"'};
%% Select limb area for analysis
l = menu('Select limb for analysis',LIMBAREAS);
%% Set directory path for finding relevant PD and Dys files and Import data
%This gets directory info for each file within the PD folder. For each file in directory, 
%need to import transcoh struct, select resting and active coh data under M1 contact 
%and place into separate (new) rest and active structures or matrices

% for PD
pdpath = uigetdir('', 'Select directory that contains PD _m1s1 files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*_transcoh.mat'); % selects '*_transcoh.mat files' that were outputs of the ecogCOH script

numPD = length(PDdir);

for i = 1:numPD
    PDcoh(i) = importdata(PDdir(i).name);
end

% for Dys
dyspath = uigetdir('', 'Select directory that contains Dys _transcoh files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*_transcoh.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numDYS = length(DYSdir);

for i = 1:numDYS
    DYScoh(i) = importdata(DYSdir(i).name);
end

% for ET
etpath = uigetdir('', 'Select directory that contains ET _transcoh files to be analyzed');
etpath = [etpath '\'];
cd(etpath);

ETdir = dir('*_transcoh.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script

numET = length(ETdir);

for i = 1:numET
    ETcoh(i) = importdata(ETdir(i).name);
end

% %% Determine which PD subjects have valid S1 contact pairs
% %If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
% for i = 1:numPD
%     S1valuesPD(i) = PDcoh(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
% end
% % A subject with contact 1 or 2 over M1 will not have an S1 contact, and
% % must be excluded from S1 analysis
% S1idx = find(S1valuesPD > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
% S1values = nonzeros(S1valuesPD); %companion to S1idx: lists the S1 contacts for those patients who have valid S1 contacts
% S1values = S1values'; 
%% Aggregate all M1, premotor, and M1-S1 coherence data from Rest state and Active state

% for PD
PD_rest = zeros(length(PDcoh(1).rest.M1S1),numPD);
PD_active = zeros(length(PDcoh(1).active.M1S1),numPD);

for i = 1: numPD
PD_rest(:,i) = PDcoh(i).rest.M1S1; %obtains rest M1S1 coherence data
PD_active(:,i) = PDcoh(i).active.M1S1; %obtains active M1S1 coherence data
end


% for Dys
DYS_rest = zeros(length(DYScoh(1).rest.M1S1),numDYS);
DYS_active = zeros(length(DYScoh(1).active.M1S1),numDYS);

for i = 1:numDYS
DYS_rest(:,i) = DYScoh(i).rest.M1S1;
DYS_active(:,i) = DYScoh(i).active.M1S1;
end

% for ET
ET_rest = zeros(length(ETcoh(1).rest.M1S1),numET);
ET_active = zeros(length(ETcoh(1).active.M1S1),numET);
 
for i = 1:numET
    ET_rest(:,i) = ETcoh(i).rest.M1S1;
    ET_active(:,i) = ETcoh(i).active.M1S1;
end

%% Average across patients
Mean_PD_rest_tcoh = mean(PD_rest,2);
Mean_PD_active_tcoh = mean(PD_active,2);
Mean_DYS_rest_tcoh = mean(DYS_rest,2);
Mean_DYS_active_tcoh = mean(DYS_active,2);
Mean_ET_rest_tcoh = mean(ET_rest,2);
Mean_ET_active_tcoh = mean(ET_active,2);

% %To plot frequency, it must be in array format, not structure format
F = PDcoh.freq; 

%% Plot by movement state
hf = figure; 
subplot(2,1,1);
plot(F, Mean_PD_rest_tcoh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_rest_tcoh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
hold on;
plot(F, Mean_ET_rest_tcoh, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
    axis([0 100 0 1]) % wide view axis
    set (gca,'xtick',[0:20:100])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')'], ['ET (n=',num2str(numET),')']);
    title ('Transformed coherence between M1-S1 during Rest')
    
subplot(2,1,2);
plot(F, Mean_PD_active_tcoh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_active_tcoh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
hold on;
plot(F, Mean_ET_active_tcoh, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
    axis([0 100 0 1]) % wide view axis
    set (gca,'xtick',[0:20:100])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')'], ['ET (n=',num2str(numET),')']);
    title (['Transformed coherence coherence between M1-S1 during' LIMBAREAS{l} 'Movement'])
%% Second Plot for M1 only
hf2 = figure;

subplot(3,1,1);
plot(F, Mean_PD_rest_tcoh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_PD_active_tcoh, 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
    axis([0 100 0 1]) % wide view axis
    set (gca,'xtick',[0:20:100])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD rest (n=',num2str(numPD),')'], ['PD active']);
    title ('Transformed coherence over M1 in PD')
    
subplot(3,1,2);
plot(F, Mean_DYS_rest_tcoh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_active_tcoh, 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
    axis([0 100 0 1]) % wide view axis
    set (gca,'xtick',[0:20:100])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['Dys rest (n=',num2str(numDYS),')'], ['Dys active']);
    title ('Transformed coherence over M1 in Dystonia')
    
subplot(3,1,3);
plot(F, Mean_ET_rest_tcoh, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_ET_active_tcoh, 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','r','LineWidth',1.5);
    axis([0 100 0 1]) % wide view axis
    set (gca,'xtick',[0:20:100])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['ET rest (n=',num2str(numET),')'], ['ET active']);
    title ('Transformed coherence over M1 in ET')
% %% Plot S1 only
% hf = figure; 
% subplot(2,1,1);
% plot(F, Mean_PD_rest_tcoh(:,:,2), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_DYS_rest_tcoh(:,:,2), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
% % hold on;
% % plot(F, Mean_ET_rest_tcoh(:,:,2), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
%     set (gca,'ytick',[0:.1:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
%     title ('Transformed coherence over S1 during Rest')
%     
% subplot(2,1,2);
% plot(F, Mean_PD_active_tcoh(:,:,2), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_DYS_active_tcoh(:,:,2), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
% % hold on;
% % plot(F, Mean_ET_active_tcoh(:,:,2), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
%     set (gca,'ytick',[0:.1:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
%     title ('Transformed coherence over S1 during Movement')
% %% Second Plot for S1 only
% hf2 = figure;
% 
% subplot(2,1,1);
% plot(F, Mean_PD_rest_tcoh(:,:,2), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_PD_active_tcoh(:,:,2), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
%     set (gca,'ytick',[0:.1:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['PD rest (n=',num2str(numPD),')'], ['PD active']);
%     title ('Transformed coherence over S1 in PD')
%     
% subplot(2,1,2);
% plot(F, Mean_DYS_rest_tcoh(:,:,2), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_DYS_active_tcoh(:,:,2), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
%     set (gca,'ytick',[0:.1:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['Dys rest (n=',num2str(numDYS),')'], ['Dys active']);
%     title ('Transformed coherence over S1 in Dystonia')
    
% subplot(3,1,3);
% plot(F, Mean_ET_rest_tcoh(:,:,2), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_ET_active_tcoh(:,:,2), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','r','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
%     set (gca,'ytick',[0:.1:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['ET rest (n=',num2str(numET),')'], ['ET active']);
%     title ('Transformed coherence over M1 in ET')
%%
% %% Plot Premotor, M1, and S1
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