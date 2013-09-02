%% timeZscore
%This program is designed to generate plots of z-scores of group timePSD and
%timeCOH data. 

%INPUT: timePSD_emg or timePSD_accel data from time_psd.m

%OUTPUT: none at this time (in future, will likely need excel file of
%z-score data)

%ALC 6/5/09

%Update: added a few lines for normalizing to baseline period. This is
%already done in the timePSD code, and in the future that code will output
%the normalized data, making this unnecessary here. 6/8/09ALC

%Update: now using code that has the data already normalized to baseline,
%so that normalization code has been commented out. 6/9/09ALC

%% Set directory path for finding relevant PD and Dys files and Import data
%This gets directory info for each file within the PD folder. For each file in directory, 
%need to import transcoh struct, select resting and active coh data under M1 contact 
%and place into separate (new) rest and active structures or matrices

% for PD
pdpath = uigetdir('', 'Select directory that contains PD _timePSD_accel or _timePSD_emg files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*.mat'); % selects all .mat files in directory **may want to update to make more specific to timePSD files

numPD = length(PDdir);

% 
for i = 1:numPD
    filename = PDdir(i).name;
    load(filename);
%     num_order = size(order, 2);
%     for j=1:num_order
    PD.dir(i) = importdata(PDdir(i).name);
end

for i = 1:numPD
    PD.dir(i).M1 = PD.dir(i).order(1);
    PD.dir(i).S1 = PD.dir(i).order(2);
    PD.dir(i).STN = PD.dir(i).order(3);
end

% Determine which PD subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numPD
    S1valuesPD(i) = PD.dir(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
end

S1idx = find(S1valuesPD > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
% S1values = isnumeric(S1valuesPD); %companion to S1idx: lists the S1 contacts for those patients who have valid S1 contacts
% S1values = S1values'; 


% noS1 = find(isnan(S1valuesPD));
% realS1valuesPD = [S1valuesPD(1:(noS1-1)) S1valuesPD((noS1+1):end)];

        % PDogram will be a large structure containing z-score analysis of
%spectrogram and coherogram data further divided by M1, S1, and STN contacts
%(represented by the digits contained in the variable "order")

% PDspect will be 3D matrix of aggregated (group) spectrogram data in the format:
%rows=freq data, col=time data, sheets = subjects
for i = 1:numPD
PD.spectrogram.M1(:,:,i) = PD.dir(i).A2plot(:,:,(PD.dir(i).M1));
PD.coherogram.M1(:,:,i) = PD.dir(i).C_trans(:,:,(PD.dir(i).M1));

PD.spectrogram.STN(:,:,i) = PD.dir(i).A2plot(:,:,(PD.dir(i).STN));
end
for i = 1:size(S1idx,2)
    PD.spectrogram.S1(:,:,i) = PD.dir(S1idx(i)).A2plot(:,:,(PD.dir(S1idx(i)).S1));
    PD.coherogram.S1(:,:,i) = PD.dir(S1idx(i)).C_trans(:,:,(PD.dir(S1idx(i)).S1));
end


% BLgrand_mean = []; % initialize structure for keeping track of mean baseline values across 
                   % patients and contacts
% for i = 1:num_chan
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    PD.spectrogram.M1mean = mean(PD.spectrogram.M1,3); 
    PD.spectrogram.M1std = std(PD.spectrogram.M1,0,3); 
    
    PD.coherogram.M1mean = mean(PD.coherogram.M1,3);
    PD.coherogram.M1std = std(PD.coherogram.M1,0,3); 
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanPD_BL = mean(PD.spectrogram.M1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    PD.spectrogram.M1bl_mean_value = mean(meanPD_BL);
    
    meanPD_BLcoh = mean(PD.coherogram.M1mean(:,1:20));
    PD.coherogram.M1bl_mean_value = mean(meanPD_BLcoh);
%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    PD.spectrogram.M1Zscore(:,:) = (PD.spectrogram.M1mean - ...
        PD.spectrogram.M1bl_mean_value) ./ PD.spectrogram.M1std;

    zcheck_psd = mean(PD.spectrogram.M1Zscore(:,1:20));
    meanBLz_psd = mean(zcheck_psd);

    PD.coherogram.M1Zscore(:,:) = (PD.coherogram.M1mean - ...
        PD.coherogram.M1bl_mean_value) ./ PD.coherogram.M1std;

    zcheck_coh = mean(PD.spectrogram.M1Zscore(:,1:20));
    meanBLz_coh = mean(zcheck_coh);

% end

% BLgrand_mean_value = mean(BLgrand_mean); % averages baseline values for all contact pairs together

% if BLgrand_mean_value ~= 1
%     disp 'The average baseline value does not equal 1!';
% end
%     for j = 1:numPD
%     zPDspect(:,:,i) = (PDogram.contact_pair(i).PDspect(:,:,j) - meanPDspect(:,:,i)) ./ ...
%         stdPDspect(:,:,i);
%     end
%% Repeat for S1
PD.spectrogram.S1mean = mean(PD.spectrogram.S1,3); 
    PD.spectrogram.S1std = std(PD.spectrogram.S1,0,3); 

PD.coherogram.S1mean = mean(PD.coherogram.S1,3); 
    PD.coherogram.S1std = std(PD.coherogram.S1,0,3); 
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanPD_BL = mean(PD.spectrogram.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    PD.spectrogram.S1bl_mean_value = mean(meanPD_BL);
    
    meanPD_BL_coh = mean(PD.coherogram.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    PD.coherogram.S1bl_mean_value = mean(meanPD_BL_coh);
%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    PD.spectrogram.S1Zscore(:,:) = (PD.spectrogram.S1mean - ...
        PD.spectrogram.S1bl_mean_value) ./ PD.spectrogram.S1std;

    zcheck_psd = mean(PD.spectrogram.S1Zscore(:,1:20));
    meanBLz_psd = mean(zcheck_psd);
    
    PD.coherogram.S1Zscore(:,:) = (PD.coherogram.S1mean - ...
        PD.coherogram.S1bl_mean_value) ./ PD.coherogram.S1std;
    
    zcheck_coh = mean(PD.coherogram.S1Zscore(:,1:20));
    meanBLz_coh = mean(zcheck_coh);
%% Repeat for STN
PD.spectrogram.STNmean = mean(PD.spectrogram.STN,3); 
    PD.spectrogram.STNstd = std(PD.spectrogram.STN,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanPD_BL = mean(PD.spectrogram.STNmean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    PD.spectrogram.STNbl_mean_value = mean(meanPD_BL);
%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    PD.spectrogram.STNZscore(:,:) = (PD.spectrogram.STNmean - ...
        PD.spectrogram.STNbl_mean_value) ./ PD.spectrogram.STNstd;

    zcheck = mean(PD.spectrogram.STNZscore(:,1:20));
    meanBLz = mean(zcheck);    
%% Repeat for DYS
dyspath = uigetdir('', 'Select directory that contains PD _timePSD_accel or _timePSD_emg files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*.mat'); % selects all .mat files in directory **may want to update to make more specific to timePSD files

numDYS = length(DYSdir);

% 
for i = 1:numDYS
    filename = DYSdir(i).name;
    load(filename);
%     num_order = size(order, 2);
%     for j=1:num_order
    DYS.dir(i) = importdata(DYSdir(i).name);
end

for i = 1:numDYS
    DYS.dir(i).M1 = DYS.dir(i).order(1);
    DYS.dir(i).S1 = DYS.dir(i).order(2);
    DYS.dir(i).STN = DYS.dir(i).order(3);
end

% Determine which DYS subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numDYS
    S1valuesDYS(i) = DYS.dir(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
end

S1idx = find(S1valuesDYS > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)

for i = 1:numDYS
DYS.spectrogram.M1(:,:,i) = DYS.dir(i).A2plot(:,:,(DYS.dir(i).M1));
DYS.coherogram.M1(:,:,i) = DYS.dir(i).C_trans(:,:,(DYS.dir(i).M1));

DYS.spectrogram.STN(:,:,i) = DYS.dir(i).A2plot(:,:,(DYS.dir(i).STN));
end
for i = 1:size(S1idx,2)
    DYS.spectrogram.S1(:,:,i) = DYS.dir(S1idx(i)).A2plot(:,:,(DYS.dir(S1idx(i)).S1));
    DYS.coherogram.S1(:,:,i) = DYS.dir(S1idx(i)).C_trans(:,:,(DYS.dir(S1idx(i)).S1));
end


% for i = 1:num_chan
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    DYS.spectrogram.M1mean = mean(DYS.spectrogram.M1,3); 
    DYS.spectrogram.M1std = std(DYS.spectrogram.M1,0,3); 
    
    DYS.coherogram.M1mean = mean(DYS.coherogram.M1,3); 
    DYS.coherogram.M1std = std(DYS.coherogram.M1,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanDYS_BL = mean(DYS.spectrogram.M1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    DYS.spectrogram.M1bl_mean_value = mean(meanDYS_BL);
    
    meanDYS_BL_coh = mean(DYS.coherogram.M1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    DYS.coherogram.M1bl_mean_value = mean(meanDYS_BL_coh);
  
    % Calculate Z score for this contact pair
    DYS.spectrogram.M1Zscore(:,:) = (DYS.spectrogram.M1mean - ...
        DYS.spectrogram.M1bl_mean_value) ./ DYS.spectrogram.M1std;

    zcheck = mean(DYS.spectrogram.M1Zscore(:,1:20));
    meanBLz = mean(zcheck);
    
    DYS.coherogram.M1Zscore(:,:) = (DYS.coherogram.M1mean - ...
        DYS.coherogram.M1bl_mean_value) ./ DYS.coherogram.M1std;

    zcheck_coh = mean(DYS.coherogram.M1Zscore(:,1:20));
    meanBLz_coh = mean(zcheck_coh);

%% Repeat for S1
DYS.spectrogram.S1mean = mean(DYS.spectrogram.S1,3); 
    DYS.spectrogram.S1std = std(DYS.spectrogram.S1,0,3); 
    
DYS.coherogram.S1mean = mean(DYS.coherogram.S1,3); 
    DYS.coherogram.S1std = std(DYS.coherogram.S1,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanDYS_BL = mean(DYS.spectrogram.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    DYS.spectrogram.S1bl_mean_value = mean(meanDYS_BL);
    
    meanDYS_BL_coh = mean(DYS.coherogram.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    DYS.coherogram.S1bl_mean_value = mean(meanDYS_BL_coh);
    
    % Calculate Z score for this contact pair
    DYS.spectrogram.S1Zscore(:,:) = (DYS.spectrogram.S1mean - ...
        DYS.spectrogram.S1bl_mean_value) ./ DYS.spectrogram.S1std;

    zcheck = mean(DYS.spectrogram.S1Zscore(:,1:20));
    meanBLz = mean(zcheck);    
    
    DYS.coherogram.S1Zscore(:,:) = (DYS.coherogram.S1mean - ...
        DYS.coherogram.S1bl_mean_value) ./ DYS.coherogram.S1std;

    zcheck_coh = mean(DYS.coherogram.S1Zscore(:,1:20));
    meanBLz_coh = mean(zcheck_coh);    

%% Repeat for STN
DYS.spectrogram.STNmean = mean(DYS.spectrogram.STN,3); 
    DYS.spectrogram.STNstd = std(DYS.spectrogram.STN,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanDYS_BL = mean(DYS.spectrogram.STNmean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    DYS.spectrogram.STNbl_mean_value = mean(meanDYS_BL);
    
    % Calculate Z score for this contact pair
    DYS.spectrogram.STNZscore(:,:) = (DYS.spectrogram.STNmean - ...
        DYS.spectrogram.STNbl_mean_value) ./ DYS.spectrogram.STNstd;

    zcheck = mean(DYS.spectrogram.STNZscore(:,1:20));
    meanBLz = mean(zcheck);    
%% Plot M1
PD2plotM1 = PD.spectrogram.M1Zscore;
DYS2plotM1 = DYS.spectrogram.M1Zscore;
 
hf1 = figure;
subplot(3,2,1);
PD2plotM1(PD2plotM1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
PD2plotM1(PD2plotM1 < (-1.96)) = (-1);
PD2plotM1(PD2plotM1~=1 & PD2plotM1~=(-1)) = 0;

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(min(PD2plotM1(:,:,:))))); % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(PD2plotM1(:,:,:)))));
clims1 = [val1 val2];

    tmp1 = PD2plotM1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf1,'textbox','String',{'Group PSD Z-scores'},'FontSize',12,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.005266 0.9639 0.1672 0.02706]);

   
    annotation(hf1,'textbox','String',{'M1'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.01874 0.8002 0.04415 0.03847]);

    annotation(hf1,'textbox','String',{'PD'},'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.2568 0.9371 0.04685 0.05436]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);

subplot(3,2,2);
DYS2plotM1(DYS2plotM1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
DYS2plotM1(DYS2plotM1 < (-1.96)) = (-1);
DYS2plotM1(DYS2plotM1~=1 & DYS2plotM1~=(-1)) = 0;

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(min(DYS2plotM1(:,:,:))))); % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(DYS2plotM1(:,:,:)))));
clims1 = [val1 val2];

    tmp1 = DYS2plotM1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf1,'textbox','String',{'DYS'},'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.7 0.9501 0.05829 0.04334]);

%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
%% Plot S1
PD2plotS1 = PD.spectrogram.S1Zscore;
DYS2plotS1 = DYS.spectrogram.S1Zscore;

subplot(3,2,3);
PD2plotS1(PD2plotS1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
PD2plotS1(PD2plotS1 < (-1.96)) = (-1);
PD2plotS1(PD2plotS1~=1 & PD2plotS1~=(-1)) = 0;
val1 = min(min(min(PD2plotS1(1:100,:,:))));
val2 = max(max(max(PD2plotS1(1:100,:,:))));
% val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
% val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = PD2plotS1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    annotation(hf1,'textbox','String',{'S1'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.01898 0.5038 0.04415 0.03847]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
    
subplot(3,2,4);
DYS2plotS1(DYS2plotS1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
DYS2plotS1(DYS2plotS1 < (-1.96)) = (-1);
DYS2plotS1(DYS2plotS1~=1 & DYS2plotS1~=(-1)) = 0;
val1 = min(min(min(DYS2plotS1(1:100,:,:))));
val2 = max(max(max(DYS2plotS1(1:100,:,:))));
% val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
% val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = DYS2plotS1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);

%% Plot STN
PD2plotSTN = PD.spectrogram.STNZscore;
DYS2plotSTN = DYS.spectrogram.STNZscore;

subplot(3,2,5);
PD2plotSTN(PD2plotSTN > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
PD2plotSTN(PD2plotSTN < (-1.96)) = (-1);
PD2plotSTN(PD2plotSTN~=1 & PD2plotSTN~=(-1)) = 0;
val1 = min(min(min(PD2plotS1(1:100,:,:))));
val2 = max(max(max(PD2plotS1(1:100,:,:))));
% val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
% val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = PD2plotSTN(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    annotation(hf1,'textbox','String',{'STN'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.01898 0.1968 0.04415 0.03847]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
    
subplot(3,2,6);
DYS2plotSTN(DYS2plotSTN > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
DYS2plotSTN(DYS2plotSTN < (-1.96)) = (-1);
DYS2plotSTN(DYS2plotSTN~=1 & DYS2plotSTN~=(-1)) = 0;
val1 = min(min(min(DYS2plotSTN(1:100,:,:))));
val2 = max(max(max(DYS2plotSTN(1:100,:,:))));
% val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
% val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = DYS2plotSTN(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');

%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);

%% Plot Coherence Z scores
% Plot M1
PD2plotM1 = PD.coherogram.M1Zscore;
DYS2plotM1 = DYS.coherogram.M1Zscore;
 
hf2 = figure;
subplot(2,2,1);
% PD2plotM1(PD2plotM1 > 1.96) = 1; 
% PD2plotM1(PD2plotM1 < (-1.96)) = (-1);
% PD2plotM1(PD2plotM1~=1 & PD2plotM1~=(-1)) = 0;

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(min(DYS2plotM1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(DYS2plotM1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = PD2plotM1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf2,'textbox','String','Group coherence Z-scores',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.003469 0.97 0.1281 0.02706]);

    annotation(hf2,'textbox','String',{'M1-STN'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.003469 0.8161 0.09267 0.04213]);
    
    annotation(hf2,'textbox','String',{'PD'},'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.2559 0.9476 0.04505 0.04335]);

%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);

subplot(2,2,2);
% DYS2plotM1(DYS2plotM1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% DYS2plotM1(DYS2plotM1 < (-1.96)) = (-1);
% DYS2plotM1(DYS2plotM1~=1 & DYS2plotM1~=(-1)) = 0;

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(min(DYS2plotM1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(DYS2plotM1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = DYS2plotM1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf2,'textbox','String',{'DYS'},'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.7135 0.952 0.04505 0.04335]);

%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);

% Plot S1
PD2plotS1 = PD.coherogram.S1Zscore;
DYS2plotS1 = DYS.coherogram.S1Zscore;

subplot(2,2,3);
% PD2plotS1(PD2plotS1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% PD2plotS1(PD2plotS1 < (-1.96)) = (-1);
% PD2plotS1(PD2plotS1~=1 & PD2plotS1~=(-1)) = 0;
% val1 = min(min(min(PD2plotS1(1:100,:,:))));
% val2 = max(max(max(PD2plotS1(1:100,:,:))));
val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = PD2plotS1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    annotation(hf2,'textbox','String',{'S1-STN'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.003703 0.3454 0.09267 0.04213]);

%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
    
subplot(2,2,4);
% DYS2plotS1(DYS2plotS1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% DYS2plotS1(DYS2plotS1 < (-1.96)) = (-1);
% DYS2plotS1(DYS2plotS1~=1 & DYS2plotS1~=(-1)) = 0;
% val1 = min(min(min(DYS2plotS1(1:100,:,:))));
% val2 = max(max(max(DYS2plotS1(1:100,:,:))));
val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = DYS2plotS1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
%     annotation(hf1, 'textbox','String','PD group coherence Z-scores from M1','HorizontalAlignment','left',...
%         'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
%% Below this is all old    
% %% Plot
% %faxis and taxis data are imported with time_psd data
% PD2plot = pdZscore;
% 
% % plot spectrogram for all ecog/lfp data
% hf1 = figure;
% % val1 = min(min(min(PD2plot(1:100,:,:))));
% % val2 = max(max(max(PD2plot(1:100,:,:))));
% val1 = (min(min(min(PD2plot(:,:,:))))) / 10; % dividing by 10 temporarily b/c min and max values so high
% val2 = (max(max(max(PD2plot(:,:,:))))) / 10;
% clims1 = [val1 val2];
% data_ch_names = {'e12','e23','e34','e45','e56','LFP'};
% 
% for i = 1:num_chan
%     subplot(2,3,i);
%     hold(gca,'on');
%     % make the time-frequency plot
%     tmp1 = PD2plot(1:100,:,i); %chopping A2plot will allow the whole colobar to be represented
%     faxis_new = faxis(1:100);
%     imagesc(taxis,faxis_new,tmp1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
% %     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
%     %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:');
%     hold(gca,'off');
%     % set the y-axis direction (YDir) to have zero at the bottom
%     set(gca,'YDir','normal');
%     % set xlim and ylim
% %     set(gca,'Xlim',[0-PRE POST]);
%     set(gca,'Ylim',[0 120]);
%     set (gca,'Xlim',[-2 2.5]);
%     % axis labels/title
%     xlabel('time (sec)');
%     ylabel('frequency (Hz)');
% %     if i==1
% % %         if typedata==1
% % % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% % %             title([outputname sprintf('\n')...
% % %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% % %                 data_ch_names{i} ' aligned to EMG Ch.' emg_ch_names{emg_i}]);
% % %         elseif typedata==2
% % % %             filename = strrep(filename,'_a.mat',''); %delete '_a.mat' from filename
% % %             title([outputname sprintf('\n')...
% % %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% % %                 data_ch_names{i} ' aligned to accel']);
% % %         else
% % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% %             title([outputname sprintf('\n')...
% %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% %                 data_ch_names{i} ' aligned to task button']);
% % %         end
% %     else
%         title(data_ch_names{i});
%         annotation(hf1, 'textbox','String','PD group PSD Z-scores','HorizontalAlignment','left',...
%         'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
% %     end
%     % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
% end
%     
% %% Repeat for Dystonia data
% % Set directory path for finding relevant Dys files and import data
% %This gets directory info for each file within the DYS folder. For each file in directory, 
% %need to import Smag_mean data (eventually A2plot data) from timePSD
% %output. Then compile and average across all patients in the group.
% 
% % for DYS
% dyspath = uigetdir('', 'Select directory that contains DYS _timePSD_accel or _timePSD_emg files to be analyzed');
% dyspath = [dyspath '\'];
% cd(dyspath);
% 
% DYSdir = dir('*.mat'); % selects all .mat files in directory **may want to update to make more specific to timePSD files
% 
% numDYS = length(DYSdir);
% 
% for i = 1:numDYS
%     filename = DYSdir(i).name;
%     load(filename);
%     
%     num_chan = size(Smag_mean, 3);
% 
%     for j=1:num_chan
% % PDogram will be a large structure containing z-score analysis of
% %spectrogram and coherogram data further divided by contact 
% % PDspect will be 3D matrix of aggregated (group) spectrogram data in the format:
% %rows=freq data, col=time data, sheets = subjects
%     
%     % ------------------------------------------------------------------
%     % temporary code to normalize each subject's data to their baseline
%     % period. In the future, this will be imported along with the other
%     % variables from _timePSD_data.mat
% %     [nfchans,nframes] = size(Smag_mean(:,:,1));
% %     first = int32((0/0.05)+1); %int32(((PRE-BL(1))/t_res)+1)
% %     last = int32(1/0.05); %int32((PRE-BL(2))/t_res)
% %         for k = 1:nfchans
% %             bl = Smag_mean(k,first:last,i);
% %             blmean = mean(bl);
% %             Smag_mean(k,:,i) = Smag_mean(k,:,i)/blmean; 
% %         end
% %     %--------------------------------------------------------------------
%     DYSogram.contact_pair(j).DYSspect(:,:,i) = A2plot(:,:,j); 
%     end
% end
% 
% for i = 1:num_chan
%     % calculate mean and std deviation for spectrogram data (power) across
%     % subjects for a given contact pair
%     DYSogram.contact_pair(i).meanDYSspect = mean(DYSogram.contact_pair(i).DYSspect,3); 
%     DYSogram.contact_pair(i).stdDYSspect = std(DYSogram.contact_pair(i).DYSspect,0,3); 
%     
%     % calculate the mean of the baseline period across subjects for given
%     % contact pair
%     meanDYS_BL = mean(DYSogram.contact_pair(i).meanDYSspect(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
%     DYSogram.contact_pair(i).BLmean_value = mean(meanDYS_BL);
% %     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
%     
%     % Want to display an error message if the baseline values are somehow
%     % not =1, but currently, this message is being displayed
%     % inappropriately (when BLmean_value = 1)
% %     if PDogram.contact_pair(i).BLmean_value ~= 1
% %        disp 'The average baseline value does not equal 1!';
% %     end
%     
%     % Calculate Z score for this contact pair
%     dysZscore(:,:,i) = (DYSogram.contact_pair(i).meanDYSspect - ...
%         DYSogram.contact_pair(i).BLmean_value) ./ DYSogram.contact_pair(i).stdDYSspect;
% 
%     zcheck = mean(dysZscore(:,1:20,i));
%     meanBLz = mean(zcheck);
% 
% end
% 
% %% Plot
% %faxis and taxis data are imported with time_psd data
% DYS2plot = dysZscore;
% 
% % plot spectrogram for all ecog/lfp data
% hf2 = figure;
% % val1 = (min(min(min(DYS2plot(1:100,:,:))))) / 10;
% % val2 = (max(max(max(DYS2plot(1:100,:,:))))) / 10;
% % clims1 = [val1 val2];
% clims1 = [-3.5 3.5]; % temporarily setting colorbar limits to these values as min/max are too big
% data_ch_names = {'e12','e23','e34','e45','e56','LFP'};
% 
% for i = 1:num_chan
%     subplot(2,3,i);
%     hold(gca,'on');
%     % make the time-frequency plot
%     tmp1 = DYS2plot(1:100,:,i); %chopping DYS2plot will allow the whole colobar to be represented
%     faxis_new = faxis(1:100);
%     imagesc(taxis,faxis_new,tmp1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
% %     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
%     %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:');
%     hold(gca,'off');
%     % set the y-axis direction (YDir) to have zero at the bottom
%     set(gca,'YDir','normal');
%     % set xlim and ylim
% %     set(gca,'Xlim',[0-PRE POST]);
%     set(gca,'Ylim',[0 120]);
%     set (gca,'Xlim',[-2 2.5]);
%     % axis labels/title
%     xlabel('time (sec)');
%     ylabel('frequency (Hz)');
% %     if i==1
% % %         if typedata==1
% % % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% % %             title([outputname sprintf('\n')...
% % %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% % %                 data_ch_names{i} ' aligned to EMG Ch.' emg_ch_names{emg_i}]);
% % %         elseif typedata==2
% % % %             filename = strrep(filename,'_a.mat',''); %delete '_a.mat' from filename
% % %             title([outputname sprintf('\n')...
% % %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% % %                 data_ch_names{i} ' aligned to accel']);
% % %         else
% % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% %             title([outputname sprintf('\n')...
% %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% %                 data_ch_names{i} ' aligned to task button']);
% % %         end
% %     else
%         title(data_ch_names{i});
%         annotation(hf2, 'textbox','String','DYS group PSD Z-scores','HorizontalAlignment','left',...
%         'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
% %     end
%     % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
% end