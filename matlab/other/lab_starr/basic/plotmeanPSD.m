% Plot mean PSD
BRAINAREAS = {'M1' 'S1' 'STN LFP'};
k = menu('Select brain area for analysis',BRAINAREAS);
%% grab PD data
% pathnamePD = uigetdir('','Select directory with PD ecogPSD.mat files');
pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numPD = length(PDdir);

restPSD = [];
activePSD=[];
counterPD = 0;

for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    if isnan(order(k));
        continue; % skip files without contact pair over selected area
    end 
    restPSD = [restPSD; rest(order(k),:)];
    activePSD = [activePSD; active(order(k),:)];

end

Yrest(:,1) = mean(restPSD);
Erest(:,1) = std(restPSD,0,1);
Yactive(:,1) = mean(activePSD);
Eactive(:,1) = std(activePSD,0,1);

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
    if isnan(order(k));
        continue; % skip files without contact pair over selected area
    end 
    DYSrest = [DYSrest; rest(order(k),:)]; %#ok<AGROW>
    DYSactive = [DYSactive; active(order(k),:)]; %#ok<AGROW>
    counterDYS = counterDYS+1;
end

Yrest(:,2) = mean(DYSrest);
Erest(:,2)= std(DYSrest,0,1);
Yactive(:,2) = mean(DYSactive);
Eactive(:,2) = std(DYSactive,0,1);

PosSD = Yrest + Erest;
NegSD = Yrest - Erest;
%% Plot line graphs with std dev
hf = figure;
plot(freq, Yrest(:,1),'color','b');
hold on;
plot(freq,Yrest(:,2),'color','g');
hold on;
plot(freq,PosSD(:,1),'color','b','linestyle','--');
hold on;
plot(freq,NegSD(:,1),'color','b','linestyle','--');
hold on;
plot(freq,PosSD(:,2),'color','g','linestyle','--');
hold on;
plot(freq,NegSD(:,2),'color','g','linestyle','--');
xlim = 100;
