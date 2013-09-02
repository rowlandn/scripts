function Ecog_GroupPSDrestabs()
% Definition: a function parallel to Ecog_GroupPSD, but for data from the 
% rest only condition which aggregates PSD data
% from each diagnostic group of interest and plots a bar graph of mean
% power in various frequency bands. The limits of these  frequency bands 
% are user defined and can be changed, however trying to plot fewer than 5 
% frequency bands can be problematic. The user also defines which brain
% area and limb is to be analyzed.

% Created by ALC 3/9/10

% INPUT = group ***_ecogPSD files from quantPSD data folder
% specifically, needs mean PSD data from the rest only condition

% CALLS function 'quantPSDrest' - this is nearly identical to the
% quantPSD subfunction called by Ecog_GroupPSD, except that the variables
% that it creates ('allfreq' and 'subfreq') have been modified to satisfy
% the rest only condition (rather than differentiating rest/active epochs).
% Calling 'quantPSDrest' can also recalcuate the allfreq and
% subfreq variables in the event that the user wishes to change the
% FREQ_QPSD variable for the purpose of analyzing group data.
% 'quantPSDrest' does not write data to an excel file.

% OUTPUT: mat file with variables for ANOVA analysis, creates plot

% 3/9/10 NOTE: this code is created in anticipation of needing group-level
% PSD data for the rest condition. As of today, the ecogPSDCOHrest code,
% which will feed into this code, does not have any output data. Thus, the
% variables herein will need to be renamed according to ecogPSDCOHrest
% outputs, once those are established. As of now, *** symbolizes places in
% the code that will likely require modifications once the ecogPSDCOHrest
% code is updated. 
%% define variables for plotting
% SAS 5/17/10: FREQ_QPSD changed to include 21.5 Hz in lower beta band
% (13-22Hz).
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

%SAS 4/28/10 Add option to run group analysis on log power
make_logplot = true;


BRAINAREAS = {'M1' 'S1' 'STN LFP'};
% LIMBAREAS = {'HAND' 'ELBOW' 'SHOULDER' 'JAW' 'FOOT' '"ARM"' '"NON-ARM"'};
% YLIM = [0 60];
%% Define frequency bands based on FREQ_QPSD 
%allows for user to create desired number of freq bins
FREQBANDS = [];
for i = 1:size(FREQ_QPSD,1)
    FREQBANDS = [FREQBANDS {num2str(FREQ_QPSD(i,:))}];
end
%% choose channel
k = menu('Select brain area for analysis',BRAINAREAS);

%% grab PD data
% pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pdpath = 'C:\Users\Starr\Documents\ECOG data\quantPSD data\PD\PSD corrected for Andrea''s study\Rest';
cd(pdpath);

PDdir = dir('*_ecogPSDrest.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numPD = length(PDdir);

PDrest = [];
counterPD = 0;

for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    if isnan(order(k));
        continue; % skip files without contact pair over selected area
    end 
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    M1_ch = order(1);
    [allfreq subfreq] = quantPSDrest(log10(psdall),FREQ_QPSD,freq,filename,M1_ch);
    
    %now create matrix of percent power in each freq band for each subject
    PDrest = [PDrest subfreq(:,1,k)]; %#ok<AGROW> %1st column contains absolute power
    counterPD = counterPD+1;
    clear psdall tcohall subfreq allfreq order 
end

Yrest(:,1) = mean(PDrest,2);
% Erest(:,1) = std(PDrest,0,2);
Erest(:,1) = std(PDrest,0,2)/(sqrt(counterPD)); %SEM

%% grab dystonia data
% dyspath = uigetdir('', 'Select directory that contains Dys _ecgPSD files to be analyzed');
dyspath = 'C:\Users\Starr\Documents\ECOG data\quantPSD data\DYS\PSD corrected for Andrea''s study\Rest';
cd(dyspath);

DYSdir = dir('*_ecogPSDrest.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
% FileListDYS = dir('*_ecogPSD.mat');

% nfilesDYT = length(FileListDYT);
numDYS = length(DYSdir);

DYSrest=[];
counterDYS=0;
for i=1:numDYS
    filename = DYSdir(i).name;
    load(filename);
    if isnan(order(k))
        continue; %skip files without contact pair over selected area
    end
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    M1_ch = order(1);
    [allfreq subfreq] = quantPSDrest(log10(psdall),FREQ_QPSD,freq,filename,M1_ch);
    
    %now create matrix of percent power in each freq band for each subject
    DYSrest = [DYSrest subfreq(:,1,k)]; %#ok<AGROW>
    counterDYS = counterDYS+1;
    clear psdall tcohall subfreq allfreq order 
end

Yrest(:,2) = mean(DYSrest,2);
% Erest(:,2)= std(DYSrest,0,2);
Erest(:,2) = std(DYSrest,0,2)/(sqrt(counterDYS)); % SEM

%% grab ET data
% if ET data needs to be analyzed, change next line to "ET=true." Otherwise
% ET=false;
ET=true;
if ET
%     etpath = uigetdir('', 'Select directory that contains ET _ecogPSD files to be analyzed');
    etpath = 'C:\Users\Starr\Documents\ECOG data\quantPSD data\ET\corrected PSD for Andrea''s study\Rest';
    cd(etpath);
    
    ETdir = dir('*_ecogPSDrest.mat'); % selects ecog PSD rest files
    numET = length(ETdir);
    
    ETrest=[];
    counterET=0;
    for i=1:numET
        filename = ETdir(i).name;
        load(filename);
        if isnan(order(k))
            continue; %skip files without contact pair over selected area
        end
         % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
    M1_ch = order(1);
    [allfreq subfreq] = quantPSDrest(log10(psdall),FREQ_QPSD,freq,filename,M1_ch);
    
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    %now create matrix of percent power in each freq band for each subject
    ETrest = [ETrest subfreq(:,1,k)]; %#ok<AGROW>
%         ETactive = [ETactive subfreq(:,4,k)]; %#ok<AGROW>
    counterET = counterET+1;
    clear psdall tcohall subfreq allfreq order   
    end

    Yrest(:,3) = mean(ETrest,2);
%     Erest(:,3) = std(ETrest,0,2);
    Erest(:,3) = std(ETrest,0,2)/(sqrt(counterET)); % SEM
end

%% grab OTHER data
%Added 4/6/09 for analysis of PD px with dystonic sx, canalso be used for
%epilepsy, etc
% if OTHER data needs to be analyzed, change next line to "OTHER=true." Otherwise
% OTHER=false;
OTHER=false;
if OTHER
%     pathnameET = uigetdir('','Select directory with OTHER ecogPSD.mat files');
    OTHERpath = uigetdir('', 'Select directory that contains OTHER _ecogPSDrest files to be analyzed');
    OTHERpath = [OTHERpath '\'];
    cd(OTHERpath);
    
    OTHERdir = dir('*_ecogPSDrest.mat'); % selects '*_ecogPSDrest.mat files'
    numOTHER = length(OTHERdir);
    
    OTHERrest=[];
    counterOTHER=0;
    for i=1:numOTHER
        filename = OTHERdir(i).name;
        load(filename);
        if isnan(order(k))
            continue; %skip files without contact pair over selected area
        end
         % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
    M1_ch = order(1);
    if exist ('P','var') % KJM epilepsy data contains a variable P which is equivalent to psdall
        psdall = P;
    end
    [allfreq subfreq] = quantPSDrest(psdall,FREQ_QPSD,freq,filename,M1_ch);
    
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    %now create matrix of percent power in each freq band for each subject
    OTHERrest = [OTHERrest subfreq(:,1,k)]; %#ok<AGROW>
    counterOTHER = counterOTHER+1;
    clear psdall tcohall subfreq allfreq order   
    end

    Yrest(:,4) = mean(OTHERrest,2);
%     Erest(:,4) = std(OTHERrest,0,2);
    Erest(:,4) = std(OTHERrest,0,2)/sqrt(counterOTHER); % SEM
end


%% plot data
hf=figure;
handlesr = barEB(Yrest,Erest);
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YMinorGrid','on');
title([BRAINAREAS{k} ' during 30 SEC REST']);
ylabel('Absolute power');
hPD=handlesr.bars(1);
hDYS=handlesr.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
if ET && ~OTHER
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['ET (n=' num2str(counterET) ')']);
    hET=handlesr.bars(3);
    set(hET,'FaceColor','m');
elseif ET && OTHER
    hET=handlesr.bars(3);
    set(hET,'FaceColor','m');
    hOTHER=handlesr.bars(4);
    set(hOTHER,'FaceColor','y');  
    
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['ET (n=' num2str(counterET) ')'], ...
        ['Epi (n=' num2str(counterOTHER) ')']);
elseif ~ET && OTHER
    hOTHER=handlesr.bars(3);
    set(hOTHER,'FaceColor','m');  
    
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['OTHER (n=' num2str(counterOTHER) ')']);
else
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')']);
%     legend(['PD'],...
%         ['Dys']);
end


%% save relevant variables
cd('C:\Users\Starr\Documents\ECOG data\quantPSD data\ANOVA stats');
outputfn=['Ecog_GroupPSDrestabs_' BRAINAREAS{k}];
save(outputfn,'PDrest','DYSrest','ETrest','Yrest','Erest');
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

% return;

