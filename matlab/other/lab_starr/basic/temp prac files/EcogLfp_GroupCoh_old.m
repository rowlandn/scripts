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
% for ET - added 1/30/09
% etpath = uigetdir('', 'Select directory that contains ET _transcoh files to be analyzed');
% etpath = [etpath '\'];
% cd(etpath);
% 
% ETdir = dir('*.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script
% 
% numET = length(ETdir);
% 
% for i = 1:numET
%     ETcoh(i) = importdata(ETdir(i).name);
% end
%% Define premotor and M1-S1 contact pairs
% PreMot_contact = M1_contact +1; %refers to contact pair anterior to M1
% S1_contact = M1_contact -2; %refers to contact pair posterior to M1
%% Aggregate all M1, premotor, and M1-S1 coherence data from Rest state and Active state

% for PD
% PD_rest = zeros(length(PDcoh(1).rest),numPD,3); % creates 3d matrix where rows=coh data, col = patients, sheets= contacts (M1, preMot, or M1S1) 
% PD_active = zeros(length(PDcoh(1).active),numPD,3);
PD_rest = zeros(length(PDcoh(1).rest),numPD,2); % creates 2D matrix for evaluating only one contact pair (M1)
PD_active = zeros(length(PDcoh(1).active),numPD,2);
for i = 1: numPD
%     PD_rest(:,i,1) = PDcoh(i).rest(:,PDcoh(i).rest_M1+1);
%     PD_rest(:,i,2) = PDcoh(i).rest(:,int8(PDcoh(i).rest_M1));
%     PD_rest(:,i,3) = PDcoh(i).rest(:,PDcoh(i).rest_M1-1);
PD_rest(:,i,1) = PDcoh(i).rest(:,int8(PDcoh(i).rest_M1)); %obtains M1 coherence data
PD_rest(:,i,2) = PDcoh(i).rest(:,int8(PDcoh(i).rest_M1) - 2); %obtains S1 coherence data
end

for i = 1:numPD
%     PD_active(:,i,1) = PDcoh(i).active(:,PDcoh(i).active_M1+1);
%     PD_active(:,i,2) = PDcoh(i).active(:,PDcoh(i).active_M1);
%     PD_active(:,i,3) = PDcoh(i).active(:,PDcoh(i).active_M1-1);
PD_active(:,i,1) = PDcoh(i).active(:,int8(PDcoh(i).active_M1));
PD_active(:,i,2) = PDcoh(i).active(:,int8(PDcoh(i).active_M1) - 2);
end

% for Dys
% DYS_rest = zeros(length(DYScoh(1).rest),numDYS,3);
% DYS_active = zeros(length(DYScoh(1).active),numDYS,3);
DYS_rest = zeros(length(DYScoh(1).rest),numDYS,2);
DYS_active = zeros(length(DYScoh(1).active),numDYS,2);
for i = 1:numDYS
%     DYS_rest(:,i,1) = DYScoh(i).rest(:,DYScoh(i).rest_M1+1);
%     DYS_rest(:,i,2) = DYScoh(i).rest(:,DYScoh(i).rest_M1);
%     DYS_rest(:,i,3) = DYScoh(i).rest(:,DYScoh(i).rest_M1-1);
DYS_rest(:,i,1) = DYScoh(i).rest(:,int8(DYScoh(i).rest_M1));
DYS_rest(:,i,2) = DYScoh(i).rest(:,int8(DYScoh(i).rest_M1) - 2);
end

for i = 1:numDYS
%     DYS_active(:,i,1) = DYScoh(i).active(:,DYScoh(i).active_M1+1);
%     DYS_active(:,i,2) = DYScoh(i).active(:,DYScoh(i).active_M1);
%     DYS_active(:,i,3) = DYScoh(i).active(:,DYScoh(i).active_M1-1);
DYS_active(:,i,1) = DYScoh(i).active(:,int8(DYScoh(i).active_M1));
DYS_active(:,i,2) = DYScoh(i).active(:,int8(DYScoh(i).active_M1) - 2);
end

% for ET
% ET_rest = zeros(length(ETcoh(1).rest),numET,2);
% ET_active = zeros(length(ETcoh(1).active),numET,2);
% 
% for i = 1:numET
%     ET_rest(:,i,1) = ETcoh(i).rest(:,int8(ETcoh(i).rest_M1));
%     ET_rest(:,i,2) = ETcoh(i).rest(:,int8(ETcoh(i).rest_M1) - 2);
% end
% 
% for i = 1:numET
%     ET_active(:,i,1) = ETcoh(i).active(:,int8(ETcoh(i).active_M1));
%     ET_active(:,i,2) = ETcoh(i).active(:,int8(ETcoh(i).active_M1) - 2);
% end
%% Average across patients
Mean_PD_rest_tcoh = mean(PD_rest,2);
Mean_PD_active_tcoh = mean(PD_active,2);
Mean_DYS_rest_tcoh = mean(DYS_rest,2);
Mean_DYS_active_tcoh = mean(DYS_active,2);
% Mean_ET_rest_tcoh = mean(ET_rest,2);
% Mean_ET_active_tcoh = mean(ET_active,2);

% %To plot frequency, it must be in array format, not structure format
F = PDcoh.freq; 

%% Plot M1 only
hf = figure; 
subplot(2,1,1);
plot(F, Mean_PD_rest_tcoh(:,:,1), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_rest_tcoh(:,:,1), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
% hold on;
% plot(F, Mean_ET_rest_tcoh(:,:,1), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
    axis([0 50 0 1]) % wide view axis
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
    title ('Transformed coherence over M1 during Rest')
    
subplot(2,1,2);
plot(F, Mean_PD_active_tcoh(:,:,1), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_active_tcoh(:,:,1), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
% hold on;
% plot(F, Mean_ET_active_tcoh(:,:,1), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
    axis([0 50 0 1]) % wide view axis
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
    title ('Transformed coherence over M1 during Movement')
%% Second Plot for M1 only
hf2 = figure;

subplot(2,1,1);
plot(F, Mean_PD_rest_tcoh(:,:,1), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_PD_active_tcoh(:,:,1), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
    axis([0 50 0 1]) % wide view axis
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD rest (n=',num2str(numPD),')'], ['PD active']);
    title ('Transformed coherence over M1 in PD')
    
subplot(2,1,2);
plot(F, Mean_DYS_rest_tcoh(:,:,1), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_active_tcoh(:,:,1), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
    axis([0 50 0 1]) % wide view axis
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['Dys rest (n=',num2str(numDYS),')'], ['Dys active']);
    title ('Transformed coherence over M1 in Dystonia')
    
% subplot(3,1,3);
% plot(F, Mean_ET_rest_tcoh(:,:,1), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','b','LineWidth',1.5);
% hold on;
% plot(F, Mean_ET_active_tcoh(:,:,1), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','r','LineWidth',1.5);
%     axis([0 50 0 1]) % wide view axis
%     set (gca,'xtick',[0:10:50])
%     set (gca,'ytick',[0:.1:1])
%     ylabel('Transformed Coherence')
%     xlabel('F')
%     legend(['ET rest (n=',num2str(numET),')'], ['ET active']);
%     title ('Transformed coherence over M1 in ET')
%% Plot S1 only
hf = figure; 
subplot(2,1,1);
plot(F, Mean_PD_rest_tcoh(:,:,2), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_rest_tcoh(:,:,2), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
% hold on;
% plot(F, Mean_ET_rest_tcoh(:,:,2), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
    axis([0 50 0 1]) % wide view axis
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
    title ('Transformed coherence over S1 during Rest')
    
subplot(2,1,2);
plot(F, Mean_PD_active_tcoh(:,:,2), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_active_tcoh(:,:,2), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','g','LineWidth',1.5);
% hold on;
% plot(F, Mean_ET_active_tcoh(:,:,2), 'DisplayName', 'ET', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in ET','color','m','LineWidth',1.5);
    axis([0 50 0 1]) % wide view axis
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD (n=',num2str(numPD),')'], ['DYS (n=' num2str(numDYS) ')']);
    title ('Transformed coherence over S1 during Movement')
%% Second Plot for S1 only
hf2 = figure;

subplot(2,1,1);
plot(F, Mean_PD_rest_tcoh(:,:,2), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_PD_active_tcoh(:,:,2), 'DisplayName', 'PD', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in PD','color','r','LineWidth',1.5);
    axis([0 50 0 1]) % wide view axis
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['PD rest (n=',num2str(numPD),')'], ['PD active']);
    title ('Transformed coherence over S1 in PD')
    
subplot(2,1,2);
plot(F, Mean_DYS_rest_tcoh(:,:,2), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','b','LineWidth',1.5);
hold on;
plot(F, Mean_DYS_active_tcoh(:,:,2), 'DisplayName', 'Dys', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence in DYS','color','r','LineWidth',1.5);
    axis([0 50 0 1]) % wide view axis
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:.1:1])
    ylabel('Transformed Coherence')
    xlabel('F')
    legend(['Dys rest (n=',num2str(numDYS),')'], ['Dys active']);
    title ('Transformed coherence over S1 in Dystonia')
    
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