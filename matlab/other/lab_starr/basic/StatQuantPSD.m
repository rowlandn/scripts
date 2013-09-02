% Calculate group stats based on allfreq/subfreq matrices

%% variables
BRAINAREAS = {'M1' 'S1' 'STN LFP'};
%% choose channel
k = menu('Select brain area for analysis',BRAINAREAS);

%% grab PD data
% pathnamePD = uigetdir('','Select directory with PD ecogPSD.mat files');
pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
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
    PDrest = [PDrest subfreq(2,2,k)];%#ok<AGROW> %2nd column contains % power at rest
    PDactive = [PDactive subfreq(2,4,k)]; %#ok<AGROW> %4th column contains % power with movement
    counterPD = counterPD+1;
end

% Yrest(:,1) = mean(PDrest,2);
% Erest(:,1) = std(PDrest,0,2);
% Yactive(:,1) = mean(PDactive,2);
% Eactive(:,1) = std(PDactive,0,2);

%% grab dystonia data
% pathnameDYT = uigetdir('','Select directory with dystonia ecogPSD.mat files');
dyspath = uigetdir('', 'Select directory that contains Dys _ecgPSD files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
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
    DYSrest = [DYSrest subfreq(2,2,k)]; %#ok<AGROW>
    DYSactive = [DYSactive subfreq(2,4,k)]; %#ok<AGROW>
    counterDYS = counterDYS+1;
end

% Yrest(:,2) = mean(DYSrest,2);
% Erest(:,2)= std(DYSrest,0,2);
% Yactive(:,2) = mean(DYSactive,2);
% Eactive(:,2) = std(DYSactive,0,2);

%% grab ET data
% if ET data needs to be analyzed, change next line to "ET=true." Otherwise
% ET=false;
ET=true;
if ET
%     pathnameET = uigetdir('','Select directory with ET ecogPSD.mat files');
    etpath = uigetdir('', 'Select directory that contains ET _ecogPSD files to be analyzed');
    etpath = [etpath '\'];
    cd(etpath);
    
    ETdir = dir('*.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script
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
        ETrest = [ETrest subfreq(2,2,k)]; %#ok<AGROW>
        ETactive = [ETactive subfreq(2,4,k)]; %#ok<AGROW>
        counterET = counterET+1;
    end

%     Yrest(:,3) = mean(ETrest,2);
%     Erest(:,3) = std(ETrest,0,2);
%     Yactive(:,3) = mean(ETactive,2);
%     Eactive(:,3) = std(ETactive,0,2);
end
%% Calculate ANOVA between dx groups
%ANOVA 1 = comparing mean percent log power in each frequency band at rest
%First create matrix containing the specific data to be compared, in the
%format rows = independent observations; columns = groups being compared
rest_beta_mtx = horzcat(PDrest', DYSrest', ETrest');
active_beta_mtx = horzcat(PDactive', DYSactive', ETactive'); 

rest_beta_mtx = 100*rest_beta_mtx;
active_beta_mtx = 100*active_beta_mtx;


[restP restTABLE] = anova1(rest_beta_mtx);
[activeP activeTABLE] = anova1(active_beta_mtx);