function Ecog_GroupPSD_abslog
% now: a function parallel to EcogLFP_GroupCOH program, which aggregates PSD data
% from each diagnostic group of interest and plots a bar graph of mean
% power in various frequency bands. The limits of these  frequency bands 
% are user defined and can be changed, however trying to plot fewer than 5 
% frequency bands can be problematic. The user also defines which brain
% area and limb is to be analyzed.

% INPUT = group _ecogPSD files from quantPSD data folder
% specifically, needs mean PSD data

% CALLS function 'quantPSDshort' - this is nearly identical to the
% quantPSD subfunction called by ecogPSD, except that it does not assign
% contact numbers or create the variable 'order' - as it already exists.
% The purpose of calling 'quantPSDshort'is to recalcuate the allfreq and
% subfreq variables in the event that the user wishes to change the
% FREQ_QPSD variable for the purpose of analyzing group data.
% 'quantPSDshort' also does not write data to an excel file.

% OUTPUT: none, creates plot

% ALC 9/16/09: changed std to sem for error bar calculation

%% define variables for plotting
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
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
% YLIM = [0 80];
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
% pdpath = [pdpath '\'];
% cd(pdpath);
pn = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pn = [pn '\'];
cd(pn);
PDdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numPD = length(PDdir);

PDrest = [];
PDactive=[];
counterPD = 0;

for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    filename = strrep(PDdir(i).name,'_ecogPSD.mat','');
    
    if isnan(order(k));
        continue; % skip files without contact pair over selected area
    end 
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename);
    [data4export] = exportdata3x31(allfreq, subfreq, filename);
    cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.
    %now create matrix of percent power in each freq band for each subject
    PDrest = [PDrest subfreq(:,1,k)];%#ok<AGROW> %2nd column contains % power at rest
    PDactive = [PDactive subfreq(:,3,k)]; %#ok<AGROW> %4th column contains % power with movement
    counterPD = counterPD+1;
end

Yrest(:,1) = mean(PDrest,2);
% Erest(:,1) = std(PDrest,0,2);
Erest(:,1) = std(PDrest,0,2)/(sqrt(counterPD));
Yactive(:,1) = mean(PDactive,2);
% Eactive(:,1) = std(PDactive,0,2);
Eactive(:,1) = std(PDactive,0,2)/(sqrt(counterPD));

%% grab dystonia data
% pathnameDYT = uigetdir('','Select directory with dystonia ecogPSD.mat files');
% dyspath = uigetdir('', 'Select directory that contains Dys _ecgPSD files to be analyzed');
% dyspath = [dyspath '\'];
% cd(dyspath);
pn = uigetdir('', 'Select directory that contains Dys _ecgPSD files to be analyzed');
pn = [pn '\'];
cd(pn);
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
    filename = strrep(DYSdir(i).name,'_ecogPSD.mat','');
    if isnan(order(k))
        continue; %skip files without contact pair over selected area
    end
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename);
    [data4export] = exportdata3x31(allfreq, subfreq, filename);
    cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.
    %now create matrix of percent power in each freq band for each subject
    DYSrest = [DYSrest subfreq(:,1,k)]; %#ok<AGROW>
    DYSactive = [DYSactive subfreq(:,3,k)]; %#ok<AGROW>
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
%     etpath = [etpath '\'];
%     cd(etpath);
    pn = uigetdir('', 'Select directory that contains ET _ecogPSD files to be analyzed');
    pn = [pn '\'];
    cd(pn); 
    ETdir = dir('*_ecogPSD.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script
    numET = length(ETdir);
    
    ETrest=[];
    ETactive=[];
    counterET=0;
    for i=1:numET
        filename = ETdir(i).name;
        load(filename);
        filename = strrep(ETdir(i).name,'_ecogPSD.mat','');
        if isnan(order(k))
            continue; %skip files without contact pair over selected area
        end
         % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
    [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename);
    [data4export] = exportdata3x31(allfreq, subfreq, filename);
    cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    %now create matrix of percent power in each freq band for each subject
    ETrest = [ETrest subfreq(:,1,k)]; %#ok<AGROW>
    ETactive = [ETactive subfreq(:,3,k)]; %#ok<AGROW>
        counterET = counterET+1;
    end

    Yrest(:,3) = mean(ETrest,2);
%     Erest(:,3) = std(ETrest,0,2);
    Erest(:,3) = std(ETrest,0,2)/(sqrt(counterET));
    Yactive(:,3) = mean(ETactive,2);
%     Eactive(:,3) = std(ETactive,0,2);
    Eactive(:,3) = std(ETactive,0,2)/(sqrt(counterET));
end

%% grab OTHER data
%Added 4/6/09 for analysis of PD px with dystonic sx
% if OTHER data needs to be analyzed, change next line to "OTHER=true." Otherwise
% ET=false;
OTHER=false;
if OTHER
%     pathnameET = uigetdir('','Select directory with OTHER ecogPSD.mat files');
    OTHERpath = uigetdir('', 'Select directory that contains OTHER _ecogPSD files to be analyzed');
    OTHERpath = [OTHERpath '\'];
    cd(OTHERpath);
    
    OTHERdir = dir('*_ecogPSD.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script
    numOTHER = length(OTHERdir);
    
    OTHERrest=[];
    OTHERactive=[];
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
    [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename);
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    %now create matrix of percent power in each freq band for each subject
    OTHERrest = [OTHERrest subfreq(:,2,k)]; %#ok<AGROW>
        OTHERactive = [OTHERactive subfreq(:,4,k)]; %#ok<AGROW>
        counterOTHER = counterOTHER+1;
    end

    Yrest(:,3) = mean(OTHERrest,2);
%     Erest(:,3) = std(OTHERrest,0,2);
    Erest(:,3) = std(OTHERrest,0,2)/sqrt(counterOTHER);
    Yactive(:,3) = mean(OTHERactive,2);
%     Eactive(:,3) = std(OTHERactive,0,2);
    Eactive(:,3) = std(OTHERactive,0,2)/sqrt(counterOTHER);
end

%% plot data
hf=figure;
%temporary code for graphing only beta band
% Yrest= Yrest(2,:);
% Erest = Erest(2,:); 
% %------
subplot(2,1,1);
handlesr = barEB(Yrest,Erest);
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YMinorGrid','on');
title([BRAINAREAS{k} ' at REST']);
ylabel('Log PSD power');
hPD=handlesr.bars(1);
hDYS=handlesr.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
if ET
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['ET (n=' num2str(counterET) ')']);
    hET=handlesr.bars(3);
    set(hET,'FaceColor','m');
elseif OTHER
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['OTHER (n=' num2str(counterOTHER) ')']);
    hOTHER=handlesr.bars(3);
    set(hOTHER,'FaceColor','m');
else
%     legend(['PD (n=' num2str(counterPD) ')'],...
%         ['Dys (n=' num2str(counterDYS) ')']);
    legend(['PD'],...
        ['Dys']);
end

subplot(2,1,2);
handlesa=barEB(Yactive,Eactive);
set(handlesa.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YMinorGrid','on');
title([BRAINAREAS{k} ' during ' LIMBAREAS{l} ' MOVEMENT']);
ylabel('Log PSD power');
hPD=handlesa.bars(1);
hDYS=handlesa.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
if ET
    hET=handlesa.bars(3);
    set(hET,'FaceColor','m');
elseif OTHER
    hOTHER=handlesa.bars(3);
    set(hOTHER,'FaceColor','m');
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

%%
function [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename)
% initialize and populate the 2x3x4 matrix 'allfreq'
% allfreq is a 2x3x4 3-D matrix in which the four 2x3 arrays correspond to the
% 4 data recording channels (premotor ecog, M1 ecog, M1-S1 ecog, and stn
% lfp). **2/19/09 - focusing on M1,S1 now, neglecting premotor for now -ALC
% The 2 rows contain:
%   1st row     -   resting state data
%   2nd row     -   active movement data
% The 3 columns contain:
%   1st col     -   frequency at which max power occurs
%   2nd col     -   max power value
%   3rd col     -   total power across all 5 frequency bands

allfreq = zeros(2,3,3);
for i=1:3
    if isnan(order(i))
        allfreq(:,:,i) = NaN;
        continue;
    end
%     rest_data = rest((order(i)),:); % from matrix-based code
%     active_data = active((order(i)),:);
    rest_data = rest.contact_pair(order(i)).log_mean_PSD;
    active_data = active.contact_pair(order(i)).log_mean_PSD;
% Need to limit rest and active data to those freq range of interest (typically 4Hz-100Hz)
%the _tmp variables are to limit the frequency spectrum for allfreq max
%values without affecting the later subfreq code, which also calls the
%rest, active, and freq variables

%to make more flexible, will use linear indexing to not require graphs to
%have 5 frequency bins
%     rest_data_tmp = rest_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
%     active_data_tmp = active_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
%     freq_tmp = freq(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));

%Need to change the _tmp such that the max freq is based on the freq bins
%and excludes, for example, 55-65Hz.
rest_data_tmp = rest_data(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));
active_data_tmp = active_data(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));
freq_tmp = freq(freq>FREQ_QPSD(1) & freq<FREQ_QPSD(end));

    [c1 i1] = max(rest_data_tmp);
    allfreq(1,1,i)=freq_tmp(i1);
    allfreq(1,2,i)=c1;
    [c2 i2] = max(active_data_tmp);
    allfreq(2,1,i)=freq_tmp(i2);
    allfreq(2,2,i)=c2;
    array1 = [];
    array2 = [];
    for j = 1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp2 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        array1 = [array1 tmp1];  %#ok<AGROW>
        array2 = [array2 tmp2]; %#ok<AGROW>
    end
    allfreq(1,3,i)=sum(array1);
    allfreq(2,3,i)=sum(array2);
end

% initialize and populate the 5x5x4 matrix 'subfreq'
% The four 5x5 arrays of subfreq correspond to the 4 data channels
% (pre-motor,M1,M1-S1,STN LFP).  Each of the 5 rows correspond to the 5
% frequency bands defined by variable FREQ_QPSD.
% The 5 columns contain:
%   1st col     -   total power in given frequency band at rest
%   2nd col     -   power in given frequency band, divided by total power
%                   in all 5 freq bands, at rest
%   3rd col     -   total power in given frequency band during movement
%   4th col     -   power in given frequency band, divided by total power
%                   in all 5 freq bands, during movement
%   5th col     -   power during movement divided by power at rest (col 3
%                   divided by col 1)
for i = 1:size(FREQ_QPSD,1)
    subfreq = zeros(i,5,3);
end

for i=1:3
    if isnan(order(i))
    subfreq(:,:,i)=NaN;
    continue;
    end
    rest_data = rest.contact_pair(order(i)).log_mean_PSD;
    active_data = active.contact_pair(order(i)).log_mean_PSD;
    for j=1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp2 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        subfreq(j,1,i) = sum(tmp1);
        subfreq(j,2,i) = sum(tmp1)/allfreq(1,3,i);
        subfreq(j,3,i) = sum(tmp2);
        subfreq(j,4,i) = sum(tmp2)/allfreq(2,3,i);
        subfreq(j,5,i) = sum(tmp2)/sum(tmp1);
    end
end
save('allfreq','subfreq');
% %% Save and write quantPSD data as excel spreadsheet
% [data4export] = exportdata3x31(allfreq, subfreq, filename);
% % cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.
% return;
function [data4export] = exportdata3x31(allfreq, subfreq, filename);
% Similar to the program "DataExport.m", the function "exportdata"
% will concatenate all the data from allfreq and
% subfreq matrices (created by quantPSD function) into a 3x31 array 
% that will be saved as an excel spreadsheet. The excel spreadsheet 
% can then be copy/pasted into SPSS. This way there is limited time 
% spent on data entry.

% 4/2/09

%% linear indexing and transposition of allfreq in preparation for concatenation 
lidxallfreq = allfreq(:);

lidxallfreq = lidxallfreq';
%% Create matrix with concatenated data from allreq and subfreq
% concatenate data into variable "data4export." This is in a format that 
% can be copy/pasted directly into the SPSS spreadsheet
data4export = [lidxallfreq(1:6) subfreq(1,:,1) subfreq(2,:,1) subfreq(3,:,1) subfreq(4,:,1) subfreq(5,:,1);...
lidxallfreq(7:12) subfreq(1,:,2) subfreq(2,:,2) subfreq(3,:,2) subfreq(4,:,2) subfreq(5,:,2);...
lidxallfreq(13:18) subfreq(1,:,3) subfreq(2,:,3) subfreq(3,:,3) subfreq(4,:,3) subfreq(5,:,3)];

%% Save data to Excel
% variable "status" tells you if it saved properly ("true")
% variable "message" will give you an error message if status=false
% If the data was not collected as Recog-Rlfp-Llimb or
% Lecog-Llfp-Rlimb, make appropriate notation in filename.

% Directory for saving the excel data 
outputpath = ['C:\Users\Starr\Documents\ECOG data analysis\Excel data\'];
cd(outputpath);

[status, message] = xlswrite([filename], data4export);

return; 
