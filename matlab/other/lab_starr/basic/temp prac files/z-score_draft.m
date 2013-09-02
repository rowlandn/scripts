%% timeZscore
%This program is designed to generate plots of z-scores of group timePSD and
%timeCOH data. 

%INPUT: timePSD_emg or timePSD_accel data from time_psd.m

%OUTPUT: none at this time (in future, will likely need excel file of
%z-score data)

%ALC 6/5/09

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

% for i = 1:numPD
%     PDogram(i) = importdata(PDdir(i).name);
% end

% % for Dys
% dyspath = uigetdir('', 'Select directory that contains Dys _timePSD_accel or _timePSD_emg files to be analyzed');
% dyspath = [dyspath '\'];
% cd(dyspath);
% 
% DYSdir = dir('*.mat'); % selects all .mat files in directory
% numDYS = length(DYSdir);
% 
% for i = 1:numDYS
%     DYSogram(i) = importdata(DYSdir(i).name);
% end
% %% Get a population mean
% % grab PD data
% pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
% pdpath = [pdpath '\'];
% cd(pdpath);
% 
% PDdir = dir('*.mat'); % selects timePSD that were outputs of the time_psd analysis script
% numPD = length(PDdir);

for i = 1:numPD
    filename = PDdir(i).name;
    load(filename);
    
    num_chan = size(Smag_mean, 3);

    for j=1:num_chan
% PDspect will be 3D matrix where rows=channels, col=PSD data, 3rd dimension = subjects        
        PDspect(j,:,i) = SMAG_mean(j,:); 
    end
end

meanPDspect = mean(PDspect,3);
stdPDspect = std(PDspect,3);
zPDspect = [];

for i=1:num_chan
    for j = 1:numPD
    zPDspect(i,:) = (PDspect(i,:,j) - meanPDspect(i,:))/stdPDspect(i,:);
    end
end

%% Plot
%faxis and taxis data are imported with time_psd data
