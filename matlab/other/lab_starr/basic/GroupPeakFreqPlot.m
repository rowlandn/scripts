%% define variables for plotting
% SAS 4/28/10: changed frequency bandwidth to be contiguous (ie.
% 4-13,13-22,22-31,31-55,76-100).  Note:21.48Hz point is included in low
% beta, which gives 5 total data points in low gamma at frequency
FREQ_QPSD = [4 13;...     % delta alpha band
             13 22;...   % low beta band
             22 31;...   % high beta band
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
YLIM = [0 80];
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
% pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pdpath = 'C:\Users\Starr\Documents\ECOG data\quantPSD data\PD\PSD corrected for Andrea''s study';
pdpath = [pdpath '\' LIMBAREAS{l}];
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

    
    %now create matrix of percent power in each freq band for each subject
    PDrest = [PDrest allfreq(1,1,k)];%#ok<AGROW> %2nd column contains % power at rest
    PDactive = [PDactive allfreq(2,1,k)]; %#ok<AGROW> %4th column contains % power with movement
    counterPD = counterPD+1;
end

Yrest(:,1) = mean(PDrest,2);
% Erest(:,1) = std(PDrest,0,2);
% Erest(:,1) = std(PDrest,0,2)/(sqrt(counterPD));
Yactive(:,1) = mean(PDactive,2);
% Eactive(:,1) = std(PDactive,0,2);
Eactive(:,1) = std(PDactive,0,2)/(sqrt(counterPD));

%% grab dystonia data
% pathnameDYT = uigetdir('','Select directory with dystonia ecogPSD.mat files');
% dyspath = uigetdir('', 'Select directory that contains Dys _ecgPSD files to be analyzed');
dyspath = 'C:\Users\Starr\Documents\ECOG data\quantPSD data\DYS\PSD corrected for Andrea''s study';
dyspath = [dyspath '\' LIMBAREAS{l}];
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
    DYSrest = [DYSrest allfreq(1,1,k)]; %#ok<AGROW>
    DYSactive = [DYSactive allfreq(2,1,k)]; %#ok<AGROW>
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
%     etpath = uigetdir('', 'Select directory that contains ET _ecogPSD files to be analyzed');
    etpath = 'C:\Users\Starr\Documents\ECOG data\quantPSD data\ET\corrected PSD for Andrea''s study';
    etpath = [etpath '\' LIMBAREAS{l}];
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

    ETrest = [ETrest allfreq(1,1,k)]; %#ok<AGROW>
    ETactive = [ETactive allfreq(2,1,k)]; %#ok<AGROW>
        counterET = counterET+1;
    end

    Yrest(:,3) = mean(ETrest,2);
%     Erest(:,3) = std(ETrest,0,2);
    Erest(:,3) = std(ETrest,0,2)/(sqrt(counterET));
    Yactive(:,3) = mean(ETactive,2);
%     Eactive(:,3) = std(ETactive,0,2);
    Eactive(:,3) = std(ETactive,0,2)/(sqrt(counterET));
end