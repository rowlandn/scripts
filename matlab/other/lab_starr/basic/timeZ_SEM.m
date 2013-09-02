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

% Minor change 3/4/10 by ALC. Now an empty structure is initiated for PD,
% DYS, or ET prior to that structure getting data put into it. Also
% eliminated code references to the variable 'order' since that variable is
% no longer used/carried through from previous analysis steps. The variable
% 'M1_ch' should now be imported and the code has been updated to use this
% variable to fill in the .M1, .S1, and .LFP sections of the structures.

% Also 3/4/10 - temporarily commenting out coherogram references because
% there is a problem with some of the coherogram data.

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
PD = [];
% 
for i = 1:numPD
    filename = PDdir(i).name;
    load(filename);
%     num_order = size(order, 2);
%     for j=1:num_order
%     PD.dir(i) = importdata(PDdir(i).name);
%       PD(i).M1 = order(1);
%       PD(i).S1 = order(2);
%       PD(i).STN = order(3);
PD(i).M1 = M1_ch;
PD(i).S1 = M1_ch - 2;
if size(A2plot,3) == 6
    PD(i).STN = 6;
else
    PD(i).STN = NaN;
end
      PD(i).A2plot = A2plot;
%       PD(i).C_trans = C_trans;
clear A2plot C_trans M1_ch
end

% for i = 1:numPD
%     PD.dir(i).M1 = PD.dir(i).order(1);
%     PD.dir(i).S1 = PD.dir(i).order(2);
%     PD.dir(i).STN = PD.dir(i).order(3);
% end

% Determine which PD subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numPD
    S1valuesPD(i) = PD(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
end

S1idx = find(S1valuesPD > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
numPDS1 = size(S1idx,2);

% Find those subjects with valid LFP data
for i = 1:numPD
    STNvaluesPD(i) = PD(i).STN;
end
STNidx = find(STNvaluesPD == 6); %since if LFP present, by definition is 6
numPDSTN = size(STNidx,2);

% PDspect will be 3D matrix of aggregated (group) spectrogram data in the format:
%rows=freq data, col=time data, sheets = subjects
for i = 1:numPD %renaming PD.spectrogram.M1 to spectrogramPD.M1, etc
spectrogramPD.M1(:,:,i) = PD(i).A2plot(:,:,(PD(i).M1));
% coherogramPD.M1(:,:,i) = PD(i).C_trans(:,:,(PD(i).M1)); - temp suspension
% spectrogramPD.STN(:,:,i) = PD(i).A2plot(:,:,(PD(i).STN)); %Not everyone has STN
end
for i = 1:size(S1idx,2) % ALC 9-10-2009: double check and make sure this is pulling the appropriate data
    spectrogramPD.S1(:,:,i) = PD(S1idx(i)).A2plot(:,:,(PD(S1idx(i)).S1));
%     coherogramPD.S1(:,:,i) = PD(S1idx(i)).C_trans(:,:,(PD(S1idx(i)).S1));
%     - temp susp
end
for i = 1:size(STNidx,2)
    spectrogramPD.STN(:,:,i) = PD(STNidx(i)).A2plot(:,:,(PD(STNidx(i)).S1));
end

% BLgrand_mean = []; % initialize structure for keeping track of mean baseline values across 
                   % patients and contacts
% for i = 1:num_chan
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    spectrogramPD.M1mean = mean(spectrogramPD.M1,3); 
    spectrogramPD.M1std = std(spectrogramPD.M1,0,3); 
    
%     coherogramPD.M1mean = mean(coherogramPD.M1,3); - temp susp
%     coherogramPD.M1std = std(coherogramPD.M1,0,3);  - temp susp

    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanPD_BL = mean(spectrogramPD.M1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    spectrogramPD.M1bl_mean_value = mean(meanPD_BL);
    
%     meanPD_BLcoh = mean(coherogramPD.M1mean(:,1:20)); - temp susp
%     coherogramPD.M1bl_mean_value = mean(meanPD_BLcoh); - temp susp

%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    spectrogramPD.M1Zscore(:,:) = (spectrogramPD.M1mean - ...
        spectrogramPD.M1bl_mean_value) ./ spectrogramPD.M1std;

    zcheck_psd = mean(spectrogramPD.M1Zscore(:,1:20));
    meanBLz_psd = mean(zcheck_psd);

%     coherogramPD.M1Zscore(:,:) = (coherogramPD.M1mean - ...
%         coherogramPD.M1bl_mean_value) ./ coherogramPD.M1std; - temp susp
% 
%     zcheck_coh = mean(coherogramPD.M1Zscore(:,1:20));
%     meanBLz_coh = mean(zcheck_coh); - temp susp

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
spectrogramPD.S1mean = mean(spectrogramPD.S1,3); 
    spectrogramPD.S1std = std(spectrogramPD.S1,0,3); 

% coherogramPD.S1mean = mean(coherogramPD.S1,3); - temp susp
%     coherogramPD.S1std = std(coherogramPD.S1,0,3);  - temp susp
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanPD_BL = mean(spectrogramPD.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    spectrogramPD.S1bl_mean_value = mean(meanPD_BL);
    
%     meanPD_BL_coh = mean(coherogramPD.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
%     coherogramPD.S1bl_mean_value = mean(meanPD_BL_coh); - temp susp of
%     these 2 lines

%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    spectrogramPD.S1Zscore(:,:) = (spectrogramPD.S1mean - ...
        spectrogramPD.S1bl_mean_value) ./ spectrogramPD.S1std;

    zcheck_psd = mean(spectrogramPD.S1Zscore(:,1:20));
    meanBLz_psd = mean(zcheck_psd);
    
%     coherogramPD.S1Zscore(:,:) = (coherogramPD.S1mean - ...
%         coherogramPD.S1bl_mean_value) ./ coherogramPD.S1std; - temp susp
%     
%     zcheck_coh = mean(coherogramPD.S1Zscore(:,1:20));
%     meanBLz_coh = mean(zcheck_coh); - temp susp
%% Repeat for STN
spectrogramPD.STNmean = mean(spectrogramPD.STN,3); 
    spectrogramPD.STNstd = std(spectrogramPD.STN,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanPD_BL = mean(spectrogramPD.STNmean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
   spectrogramPD.STNbl_mean_value = mean(meanPD_BL);
%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    spectrogramPD.STNZscore(:,:) = (spectrogramPD.STNmean - ...
        spectrogramPD.STNbl_mean_value) ./ spectrogramPD.STNstd;

    zcheck = mean(spectrogramPD.STNZscore(:,1:20));
    meanBLz = mean(zcheck);    
%% Repeat for DYS
dyspath = uigetdir('', 'Select directory that contains DYS _timePSD_accel or _timePSD_emg files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*.mat'); % selects all .mat files in directory **may want to update to make more specific to timePSD files

numDYS = length(DYSdir);
DYS = [];
% 
for i = 1:numDYS
    filename = DYSdir(i).name;
    load(filename);
%     num_order = size(order, 2);
%     for j=1:num_order
%     DYS.dir(i) = importdata(DYSdir(i).name);
%       DYS(i).M1 = order(1);
%       DYS(i).S1 = order(2);
%       DYS(i).STN = order(3);
DYS(i).M1 = M1_ch;
DYS(i).S1 = M1_ch - 2;
if size(A2plot,3) == 6
    DYS(i).STN = 6;
else
    DYS(i).STN = NaN;
end
      DYS(i).A2plot = A2plot;
%       DYS(i).C_trans = C_trans;
clear A2plot C_trans M1_ch
end

% for i = 1:numDYS
%     DYS.dir(i).M1 = DYS.dir(i).order(1);
%     DYS.dir(i).S1 = DYS.dir(i).order(2);
%     DYS.dir(i).STN = DYS.dir(i).order(3);
% end

% Determine which DYS subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numDYS
    S1valuesDYS(i) = DYS(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
end

S1idx = find(S1valuesDYS > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
numDYSS1 = size(S1idx,2);

% Find those subjects with valid LFP data
for i = 1:numDYS
    STNvaluesDYS(i) = DYS(i).STN;
end
STNidx = find(STNvaluesDYS == 6); %since if LFP present, by definition is 6
numDYSSTN = size(STNidx,2);

for i = 1:numDYS
spectrogramDYS.M1(:,:,i) = DYS(i).A2plot(:,:,(DYS(i).M1));
% coherogramDYS.M1(:,:,i) = DYS(i).C_trans(:,:,(DYS(i).M1)); - temp susp

% DYS.spectrogram.STN(:,:,i) = DYS.dir(i).A2plot(:,:,(DYS.dir(i).STN));
end
for i = 1:size(S1idx,2)
    spectrogramDYS.S1(:,:,i) = DYS(S1idx(i)).A2plot(:,:,(DYS(S1idx(i)).S1));
%     coherogramDYS.S1(:,:,i) = DYS(S1idx(i)).C_trans(:,:,(DYS(S1idx(i)).S1)); - temp susp 
end
for i = 1:size(STNidx,2)
    spectrogramDYS.STN(:,:,i) = DYS(STNidx(i)).A2plot(:,:,(DYS(STNidx(i)).S1));
end

% for i = 1:num_chan
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    spectrogramDYS.M1mean = mean(spectrogramDYS.M1,3); 
    spectrogramDYS.M1std = std(spectrogramDYS.M1,0,3); 
    
%     coherogramDYS.M1mean = mean(coherogramDYS.M1,3); - temp susp
%     coherogramDYS.M1std = std(coherogramDYS.M1,0,3); - temp susp
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanDYS_BL = mean(spectrogramDYS.M1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    spectrogramDYS.M1bl_mean_value = mean(meanDYS_BL);
    
%     meanDYS_BL_coh = mean(coherogramDYS.M1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
%     coherogramDYS.M1bl_mean_value = mean(meanDYS_BL_coh); - these 2 lines
%     temp susp
  
    % Calculate Z score for this contact pair
    spectrogramDYS.M1Zscore(:,:) = (spectrogramDYS.M1mean - ...
        spectrogramDYS.M1bl_mean_value) ./ spectrogramDYS.M1std;

    zcheck = mean(spectrogramDYS.M1Zscore(:,1:20));
    meanBLz = mean(zcheck);
    
%     coherogramDYS.M1Zscore(:,:) = (coherogramDYS.M1mean - ...
%         coherogramDYS.M1bl_mean_value) ./ coherogramDYS.M1std; -temp susp
% 
%     zcheck_coh = mean(coherogramDYS.M1Zscore(:,1:20));
%     meanBLz_coh = mean(zcheck_coh); - temp susp

%% Repeat for S1
spectrogramDYS.S1mean = mean(spectrogramDYS.S1,3); 
    spectrogramDYS.S1std = std(spectrogramDYS.S1,0,3); 
    
% coherogramDYS.S1mean = mean(coherogramDYS.S1,3); - temp susp
%     coherogramDYS.S1std = std(coherogramDYS.S1,0,3);  - temp susp
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanDYS_BL = mean(spectrogramDYS.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    spectrogramDYS.S1bl_mean_value = mean(meanDYS_BL);
    
%     meanDYS_BL_coh = mean(coherogramDYS.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
%     coherogramDYS.S1bl_mean_value = mean(meanDYS_BL_coh); - these 2 lines
%     temp susp
    
    % Calculate Z score for this contact pair
    spectrogramDYS.S1Zscore(:,:) = (spectrogramDYS.S1mean - ...
        spectrogramDYS.S1bl_mean_value) ./ spectrogramDYS.S1std;

    zcheck = mean(spectrogramDYS.S1Zscore(:,1:20));
    meanBLz = mean(zcheck);    
    
%     coherogramDYS.S1Zscore(:,:) = (coherogramDYS.S1mean - ...
%         coherogramDYS.S1bl_mean_value) ./ coherogramDYS.S1std; -temp susp
% 
%     zcheck_coh = mean(coherogramDYS.S1Zscore(:,1:20));
%     meanBLz_coh = mean(zcheck_coh);    - temp susp

%% Repeat for STN
spectrogramDYS.STNmean = mean(spectrogramDYS.STN,3); 
    spectrogramDYS.STNstd = std(spectrogramDYS.STN,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanDYS_BL = mean(spectrogramDYS.STNmean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    spectrogramDYS.STNbl_mean_value = mean(meanDYS_BL);
    
    % Calculate Z score for this contact pair
    spectrogramDYS.STNZscore(:,:) = (spectrogramDYS.STNmean - ...
        spectrogramDYS.STNbl_mean_value) ./ spectrogramDYS.STNstd;

    zcheck = mean(spectrogramDYS.STNZscore(:,1:20));
    meanBLz = mean(zcheck);    
    
%% Repeat for ET
% ET = true;
%   ET = false;
etpath = uigetdir('', 'Select directory that contains ET _timePSD_accel or _timePSD_emg files to be analyzed');
etpath = [etpath '\'];
cd(etpath);

ETdir = dir('*.mat'); % selects all .mat files in directory **may want to update to make more specific to timePSD files

numET = length(ETdir);
ET = [];
% 
for i = 1:numET
    filename = ETdir(i).name;
    load(filename);
%     num_order = size(order, 2);
%     for j=1:num_order
%     DYS.dir(i) = importdata(DYSdir(i).name);
%       ET(i).M1 = order(1);
%       ET(i).S1 = order(2);
%       ET(i).STN = order(3);
      ET(i).M1 = M1_ch;
      ET(i).S1 = M1_ch - 2; 
      ET(i).A2plot = A2plot;
%       ET(i).C_trans = C_trans;
clear A2plot C_trans M1_ch
end

% for i = 1:numDYS
%     DYS.dir(i).M1 = DYS.dir(i).order(1);
%     DYS.dir(i).S1 = DYS.dir(i).order(2);
%     DYS.dir(i).STN = DYS.dir(i).order(3);
% end

% Determine which DYS subjects have valid S1 contact pairs
%If M1 contact is contact 1 or 2, the patient does not have a valid S1 pair
for i = 1:numET
    S1valuesET(i) = ET(i).S1; %defines S1 contact for each subject based on that subject's M1 contact
end

S1idx = find(S1valuesET > 0); %finds positive, non-zero elements of S1values (lists subjects with valid S1 contacts)
numETS1 = size(S1idx,2);

for i = 1:numET
spectrogramET.M1(:,:,i) = ET(i).A2plot(:,:,(ET(i).M1));
% coherogramDYS.M1(:,:,i) = DYS(i).C_trans(:,:,(DYS(i).M1));

% DYS.spectrogram.STN(:,:,i) = DYS.dir(i).A2plot(:,:,(DYS.dir(i).STN));
end
for i = 1:size(S1idx,2)
    spectrogramET.S1(:,:,i) = ET(S1idx(i)).A2plot(:,:,(ET(S1idx(i)).S1));
%     coherogramDYS.S1(:,:,i) = DYS(S1idx(i)).C_trans(:,:,(DYS(S1idx(i)).S1));
end
% for i = 1:size(STNidx,2)
%     spectrogramDYS.STN(:,:,i) = DYS(STNidx(i)).A2plot(:,:,(DYS(STNidx(i)).S1));
% end

% for i = 1:num_chan
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    spectrogramET.M1mean = mean(spectrogramET.M1,3); 
    spectrogramET.M1std = std(spectrogramET.M1,0,3); 
    
%     coherogramDYS.M1mean = mean(coherogramDYS.M1,3); 
%     coherogramDYS.M1std = std(coherogramDYS.M1,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanET_BL = mean(spectrogramET.M1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    spectrogramET.M1bl_mean_value = mean(meanET_BL);
    
%     meanDYS_BL_coh = mean(coherogramDYS.M1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
%     coherogramDYS.M1bl_mean_value = mean(meanDYS_BL_coh);
  
    % Calculate Z score for this contact pair
    spectrogramET.M1Zscore(:,:) = (spectrogramET.M1mean - ...
        spectrogramET.M1bl_mean_value) ./ spectrogramET.M1std;

    zcheck = mean(spectrogramET.M1Zscore(:,1:20));
    meanBLz = mean(zcheck);
    
%     coherogramDYS.M1Zscore(:,:) = (coherogramDYS.M1mean - ...
%         coherogramDYS.M1bl_mean_value) ./ coherogramDYS.M1std;

%     zcheck_coh = mean(coherogramDYS.M1Zscore(:,1:20));
%     meanBLz_coh = mean(zcheck_coh);

%% Repeat for S1
spectrogramET.S1mean = mean(spectrogramET.S1,3); 
    spectrogramET.S1std = std(spectrogramET.S1,0,3); 
    
% coherogramET.S1mean = mean(coherogramET.S1,3); 
%     coherogramET.S1std = std(coherogramET.S1,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanET_BL = mean(spectrogramET.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    spectrogramET.S1bl_mean_value = mean(meanET_BL);
    
%     meanDYS_BL_coh = mean(coherogramDYS.S1mean(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
%     coherogramDYS.S1bl_mean_value = mean(meanDYS_BL_coh);
    
    % Calculate Z score for this contact pair
    spectrogramET.S1Zscore(:,:) = (spectrogramET.S1mean - ...
        spectrogramET.S1bl_mean_value) ./ spectrogramET.S1std;

    zcheck = mean(spectrogramET.S1Zscore(:,1:20));
    meanBLz = mean(zcheck);    
    
%     coherogramDYS.S1Zscore(:,:) = (coherogramDYS.S1mean - ...
%         coherogramDYS.S1bl_mean_value) ./ coherogramDYS.S1std;

%     zcheck_coh = mean(coherogramDYS.S1Zscore(:,1:20));
%     meanBLz_coh = mean(zcheck_coh);  
%% Plot M1
% Display spectrogram data for relevant frequencies only (1-100Hz)
PD2plotM1 = spectrogramPD.M1Zscore(1:120,:); 
DYS2plotM1 = spectrogramDYS.M1Zscore(1:120,:);
ET2plotM1 = spectrogramET.M1Zscore(1:120,:);
 
hf1 = figure;
subplot(3,3,1);

% Lump values by std dev 
PD2plotM1(PD2plotM1 >= 2) = 2;
PD2plotM1(1<=PD2plotM1 & PD2plotM1<2) = 1;
PD2plotM1(-1< PD2plotM1 & PD2plotM1 < 1) = 0;
PD2plotM1((-2) < PD2plotM1 & PD2plotM1 <= (-1)) = (-1);
PD2plotM1(PD2plotM1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
% val1 = (min(min(PD2plotM1(:,:,:)))); 
% val2 = (max(max(PD2plotM1(:,:,:))));
% clims1 = [val1 val2];
clims1 = [-2 2];
    imagesc(taxis,faxis,PD2plotM1,clims1); 
    colormap(gray); %makes plot black and white (greyscale)
  
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % define tick marks along y-axis
    set(gca,'YTick',(0:20:120));


hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:'); 

    % axis labels/title  
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    % list n for each group/brain area plot
%     annotation(hf1,'textbox','String',{'n=' numPD},'FontName',...
%         'Arial','FontSize',9,'FitBoxToText','off','LineStyle','none',...
%         'Position',[0.299 0.927 0.041 0.023]);
        
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
    'Position',[0.211 0.948 0.058 0.043]);
%     end
    % put a color scale indicator next to the time-coherence plot
    colorbar([0.9307 0.1048 0.02354 0.8226]);

subplot(3,3,2);
% DYS2plotM1(DYS2plotM1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% DYS2plotM1(DYS2plotM1 < (-1.96)) = (-1);
% DYS2plotM1(DYS2plotM1~=1 & DYS2plotM1~=(-1)) = 0;

DYS2plotM1(DYS2plotM1 >= 2) = 2;
DYS2plotM1(1<=DYS2plotM1 & DYS2plotM1<2) = 1;
DYS2plotM1(-1< DYS2plotM1 & DYS2plotM1 < 1) = 0;
DYS2plotM1((-2) < DYS2plotM1 & DYS2plotM1 <= (-1)) = (-1);
DYS2plotM1(DYS2plotM1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

val1 = (min(min(DYS2plotM1(:,:,:)))); 
val2 = (max(max(DYS2plotM1(:,:,:))));
% clims1 = [val1 val2];

% Plot data
imagesc(taxis,faxis,DYS2plotM1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    set(gca,'YTick',(0:20:120));
    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    % list n for each group/brain area plot
%     annotation(hf1,'textbox','String',{'n =' num2str(numDYS)},'FontName',...
%         'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
%         'Position',[0.579 0.927 0.041 0.023]);
    
    annotation(hf1,'textbox','String',{'DYS'},'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.488  0.948 0.058 0.043]);

subplot(3,3,3);
% ET2plotM1(ET2plotM1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% ET2plotM1(ET2plotM1 < (-1.96)) = (-1);
% ET2plotM1(ET2plotM1~=1 & ET2plotM1~=(-1)) = 0;

ET2plotM1(ET2plotM1 >= 2) = 2;
ET2plotM1(1<=ET2plotM1 & ET2plotM1<2) = 1;
ET2plotM1(-1< ET2plotM1 & ET2plotM1 < 1) = 0;
ET2plotM1((-2) < ET2plotM1 & ET2plotM1 <= (-1)) = (-1);
ET2plotM1(ET2plotM1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(ET2plotM1(:,:,:)))); % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(ET2plotM1(:,:,:))));
% clims1 = [val1 val2];

imagesc(taxis,faxis,ET2plotM1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % define Y-axis tickmarks
    set(gca,'YTick',[0:20:120]);

    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    % list n for each group/brain area plot   
%     annotation(hf1,'textbox','String',{'n =' num2str(numET)},'FontName',...
%         'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
%         'Position',[0.859 0.927 0.041 0.023]);
    
    annotation(hf1,'textbox','String',{'ET'},'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.771  0.948 0.058 0.043]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
%% Plot S1
PD2plotS1 = spectrogramPD.S1Zscore(1:120,:);
DYS2plotS1 = spectrogramDYS.S1Zscore(1:120,:);
ET2plotS1 = spectrogramET.S1Zscore(1:120,:);

subplot(3,3,4);
PD2plotS1(PD2plotS1 >= 2) = 2;
PD2plotS1(1<=PD2plotS1 & PD2plotS1<2) = 1;
PD2plotS1(-1< PD2plotS1 & PD2plotS1 < 1) = 0;
PD2plotS1((-2) < PD2plotS1 & PD2plotS1 <= (-1)) = (-1);
PD2plotS1(PD2plotS1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = min(min(min(PD2plotS1(1:100,:,:))));
% val2 = max(max(max(PD2plotS1(1:100,:,:))));
% % val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
% % val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
% clims1 = [val1 val2];
clims1 = [-2 2];

imagesc(taxis,faxis,PD2plotS1,clims1);
    colormap(gray); %makes plot black and white (greyscale)

    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    set(gca,'YTick',(0:20:120));
    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    annotation(hf1,'textbox','String',{'S1'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.01898 0.5038 0.04415 0.03847]);
    % list n for each group/brain area plot
%     annotation(hf1,'textbox','String',{'n =' num2str(numPDS1)},'FontName',...
%         'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
%         'Position',[0.299 0.628 0.041 0.023]);

    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
    
subplot(3,3,5);
DYS2plotS1(DYS2plotS1 >= 2) = 2;
DYS2plotS1(1<=DYS2plotS1 & DYS2plotS1<2) = 1;
DYS2plotS1(-1< DYS2plotS1 & DYS2plotS1 < 1) = 0;
DYS2plotS1((-2) < DYS2plotS1 & DYS2plotS1 <= (-1)) = (-1);
DYS2plotS1(DYS2plotS1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = (min(min(DYS2plotS1(:,:,:)))); 
% val2 = (max(max(DYS2plotS1(:,:,:))));
% clims1 = [val1 val2];

% Plot data
imagesc(taxis,faxis,DYS2plotS1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    set(gca,'YTick',(0:20:120));
    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    % list n for each group/brain area plot
%     annotation(hf1,'textbox','String',{'n =' num2str(numDYSS1)},'FontName',...
%         'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
%         'Position',[0.579 0.628 0.041 0.023]);
    
subplot(3,3,6);
ET2plotS1(ET2plotS1 >= 2) = 2;
ET2plotS1(1<=ET2plotS1 & ET2plotS1<2) = 1;
ET2plotS1(-1< ET2plotS1 & ET2plotS1 < 1) = 0;
ET2plotS1((-2) < ET2plotS1 & ET2plotS1 <= (-1)) = (-1);
ET2plotS1(ET2plotS1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];


% val1 = (min(min(ET2plotS1(:,:,:)))); 
% val2 = (max(max(ET2plotS1(:,:,:))));
% clims1 = [val1 val2];

imagesc(taxis,faxis,ET2plotS1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % define Y-axis tickmarks
    set(gca,'YTick',[0:20:120]);

    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    % list n for each group/brain area plot    
%     annotation(hf1,'textbox','String',{'n =' num2str(numETS1)},'FontName',...
%         'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
%         'Position',[0.859 0.628 0.041 0.023]);
%% Plot STN
PD2plotSTN = spectrogramPD.STNZscore(1:120,:);
DYS2plotSTN = spectrogramDYS.STNZscore(1:120,:);

subplot(3,3,7);
PD2plotSTN(PD2plotSTN >= 2) = 2;
PD2plotSTN(1<=PD2plotSTN & PD2plotSTN<2) = 1;
PD2plotSTN(-1< PD2plotSTN & PD2plotSTN < 1) = 0;
PD2plotSTN((-2) < PD2plotSTN & PD2plotSTN <= (-1)) = (-1);
PD2plotSTN(PD2plotSTN <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
% val1 = (min(min(PD2plotM1(:,:,:)))); 
% val2 = (max(max(PD2plotM1(:,:,:))));
% clims1 = [val1 val2];
clims1 = [-2 2];
    imagesc(taxis,faxis,PD2plotSTN,clims1); 
    colormap(gray); %makes plot black and white (greyscale)
  
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % define tick marks along y-axis
    set(gca,'YTick',(0:20:120));


hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:'); 

    % axis labels/title  
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    % list n for each group/brain area plot
%     annotation(hf1,'textbox','String',{'n =' num2str(numPDSTN)},'FontName',...
%         'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
%         'Position',[0.299 0.328 0.041 0.023]);
    
    annotation(hf1,'textbox','String',{'STN'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.01898 0.1968 0.04415 0.03847]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
    
subplot(3,3,8);
DYS2plotSTN(DYS2plotSTN >= 2) = 2;
DYS2plotSTN(1<=DYS2plotSTN & DYS2plotSTN<2) = 1;
DYS2plotSTN(-1< DYS2plotSTN & DYS2plotSTN < 1) = 0;
DYS2plotSTN((-2) < DYS2plotSTN & DYS2plotSTN <= (-1)) = (-1);
DYS2plotSTN(DYS2plotSTN <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = (min(min(DYS2plotSTN(:,:,:)))); 
% val2 = (max(max(DYS2plotSTN(:,:,:))));
% clims1 = [val1 val2];

% Plot data
imagesc(taxis,faxis,DYS2plotSTN,clims1);
    colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    set(gca,'YTick',(0:20:120));
    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    % list n for each group/brain area plot
%     annotation(hf1,'textbox','String',{'n =' num2str(numDYSSTN)},'FontName',...
%         'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
%         'Position',[0.579 0.328 0.041 0.023]);
    
% subplot(3,3,9);
    
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);