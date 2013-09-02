function PlotQuantPSD
% plot results of quantPSD output from ecogPSD analysis

%%define variables
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   
FREQBANDS = {'4-12' '13-21' '22-30' '31-55' '76-100'};
AREAS = {'pre-motor' 'M1' 'S1' 'STN LFP'};
YLIM = [0 70];
%% choose channel
k = menu('Select area for analysis',AREAS);

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
    DYSrest = [DYSrest subfreq(:,2,k)]; %#ok<AGROW>
    DYSactive = [DYSactive subfreq(:,4,k)]; %#ok<AGROW>
    counterDYS = counterDYS+1;
end

Yrest(:,1) = mean(DYSrest,2);
Erest(:,1)= std(DYSrest,0,2);
Yactive(:,1) = mean(DYSactive,2);
Eactive(:,1) = std(DYSactive,0,2);

%% grab PD data
% pathnamePD = uigetdir('','Select directory with PD ecogPSD.mat files');
% pathnamePD='C:\Documents and Settings\Sho\My Documents\Lab Documents\Ecog-LFP-SSEP\DMRF grant\PD elbow';
% cd(pathnamePD);
% FileListPD = dir('*_ecogPSD.mat');
% nfilesPD = length(FileListPD);
PDpath = uigetdir('', 'Select directory that contains PD _ecgPSD files to be analyzed');
PDpath = [PDpath '\'];
cd(PDpath);

PDdir = dir('*.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
% FileListDYS = dir('*_ecogPSD.mat');

% nfilesDYT = length(FileListDYT);
numPD = length(DYSdir);
PDrest = [];
PDactive=[];
counterPD = 0;
for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    if isnan(order(k));
        continue; % skip files without contact pair over selected area
    end
    PDrest = [PDrest subfreq(:,2,k)];%#ok<AGROW> %2nd column contains % power at rest
    PDactive = [PDactive subfreq(:,4,k)]; %#ok<AGROW> %4th column contains % power with movement
    counterPD = counterPD+1;
end

Yrest(:,2) = mean(PDrest,2);
Erest(:,2) = std(PDrest,0,2);
Yactive(:,2) = mean(PDactive,2);
Eactive(:,2) = std(PDactive,0,2);

%% grab ET data
% if ET data needs to be analyzed, change next line to "ET=true." Otherwise
% "ET=false"
ET=false;
if ET
%     pathnameET = uigetdir('','Select directory with ET ecogPSD.mat files');
    pathnameET = 'C:\Documents and Settings\Sho\My Documents\Lab Documents\Ecog-LFP-SSEP\DMRF grant\ET elbow';
    cd(pathnameET);
    FileListET = dir('*_ecogPSD.mat');
    nfilesET = length(FileListET);
    ETrest=[];
    ETactive=[];
    counterET=0;
    for i=1:nfilesET
        filename = FileListET(i).name;
        load(filename);
        if isnan(order(k))
            continue; %skip files without contact pair over selected area
        end
        ETrest = [ETrest subfreq(:,2,k)]; %#ok<AGROW>
        ETactive = [ETactive subfreq(:,4,k)]; %#ok<AGROW>
        counterET = counterET+1;
    end

    Yrest(:,3) = mean(ETrest,2);
    Erest(:,3) = std(ETrest,0,2);
    Yactive(:,3) = mean(ETactive,2);
    Eactive(:,3) = std(ETactive,0,2);
end
%% plot data
hf=figure;

subplot(2,1,1);
handlesr = barEB(100*Yrest,100*Erest);
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YMinorGrid','on');
title([AREAS{k} ' at REST']);
ylabel('% power');
hPD=handlesr.bars(1);
hDYT=handlesr.bars(2);
set(hPD,'FaceColor','b');
set(hDYT,'FaceColor','r');
if ET
    legend(['Dyst (n=' num2str(counterDYT) ')'],...
        ['PD (n=' num2str(counterPD) ')'],...
        ['ET (n=' num2str(counterET) ')']);
    hET=handlesr.bars(3);
    set(hET,'FaceColor','y');
else
    legend(['Dyst (n=' num2str(counterDYS) ')'],...
        ['PD (n=' num2str(counterPD) ')']);
end

subplot(2,1,2);
handlesa=barEB(100*Yactive,100*Eactive);
set(handlesa.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YMinorGrid','on');
title([AREAS{k} ' during MOVEMENT']);
ylabel('% power');
hPD=handlesa.bars(1);
hDYT=handlesa.bars(2);
set(hPD,'FaceColor','b');
set(hDYT,'FaceColor','r');
if ET
    hET=handlesa.bars(3);
    set(hET,'FaceColor','y');
end
function handles = barEB(barvals,errorvals)

% barEB plots bar graph of m-by-n matrix of bar values (by calling MATLAB
% function "bar") as well as error bars (by calling "errorbar")

% note: barvals and errorvals must be the same size!

% a list of handles are listed so that the user can change the properties
% of the plot
% handles.bars: handle to bar plot
% handles.gca: handle to axis

[numgroups numbars] = size(barvals); % # of groups, # bars in each group

handles.bars = bar(barvals);
hold on;
groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    handles.errors(i) = errorbar(x, barvals(:,i), errorvals(:,i), 'k', 'linestyle', 'none');
end

handles.ca = gca;
hold off;
return;