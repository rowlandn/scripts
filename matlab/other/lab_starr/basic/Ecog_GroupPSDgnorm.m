function Ecog_GroupPSDgnorm
% Definition: same as Ecog_GroupPSD, except bars represent the log power
% in each frequency band normalized by high gamma power, rather than the percent power relative to the
% total power across the five frequency bands - ALC 3/9/10

% INPUT = group _ecogPSD files from quantPSD data folder
% specifically, needs log mean PSD data

% CALLS function 'quantPSDshort_loggnorm' - this is nearly identical to the
% quantPSD subfunction called by ecogPSD, except that it does not assign
% contact numbers or create the variable 'order' - as it already exists.
% The purpose of calling 'quantPSDshort_loggnorm'is to recalcuate the allfreq and
% subfreq variables in the event that the user wishes to change the
% FREQ_QPSD variable for the purpose of analyzing group data.
% 'quantPSDshort_loggnorm' also does not write data to an excel file.

% OUTPUT: mat file with variables for ANOVA analysis, creates 2 plots, one for low/high beta and another for each
% frequency band
%% define variables for plotting
% SAS 4/29/10: changed frequency bandwidth so that they are contiguous (ie.
% 4-13,13-22,22-31,31-55,76-100).  Note:21.48Hz point is included in low
% beta, which gives 5 total data points in low gamma at frequency
% resolution of 1.95 Hz.
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
pdpath = 'C:\Users\Starr\Documents\ECOG data\quantPSD data\PD\PSD corrected for Andrea''s study';
pdpath = [pdpath '\' LIMBAREAS{l}];
cd(pdpath);

PDdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numPD = length(PDdir);

% PDrest and PDactive will contain subfrequency power
PDrest = [];
PDactive = [];
% PDrestbeta and PDactivebeta will contain just beta power
PDrestbeta = [];
PDactivebeta=[];

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
    %     [allfreq subfreq] =
    %     quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);

    [norm_rest_beta norm_active_beta subfreq] = ...
        quantPSDshort_loggnorm(rest,active,FREQ_QPSD,order,freq,filename);
    
    %now create matrix of percent power in each freq band for each subject
    PDrest = [PDrest subfreq(:,2,k)];%#ok<AGROW> %2nd column contains power at rest normalized by gamma power at rest
    PDactive = [PDactive subfreq(:,4,k)]; %#ok<AGROW> %4th column contains power with movement normalized by gamma power with movement
    
    PDrestbeta = [PDrestbeta norm_rest_beta(k)]; %#ok<AGROW> 
    PDactivebeta = [PDactivebeta norm_active_beta(k)]; %#ok<AGROW> 
    
    counterPD = counterPD+1;
end

Yrest(:,1) = mean(PDrest,2);
Erest(:,1) = std(PDrest,0,2)/(sqrt(counterPD));
Yactive(:,1) = mean(PDactive,2);
Eactive(:,1) = std(PDactive,0,2)/(sqrt(counterPD));

% Ybeta is a 2x3 matrix.  2 rows for rest and active.  3 columns for PD, DYS,
% and ET
Ybeta(1,1) = mean(PDrestbeta);
Ebeta(1,1) = std(PDrestbeta)/(sqrt(counterPD));
Ybeta(2,1) = mean(PDactivebeta);
Ebeta(2,1) = std(PDactivebeta)/(sqrt(counterPD));

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

%DYSrest and DYSactive will contain subfrequency power
DYSrest = [];
DYSactive = [];
%DYSrestbeta and DYSactivebeta will contain beta power
DYSrestbeta=[];
DYSactivebeta=[];

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
    
    [norm_rest_beta norm_active_beta subfreq] = ...
        quantPSDshort_loggnorm(rest,active,FREQ_QPSD,order,freq,filename);
    
    %now create matrix of percent power in each freq band for each subject
    DYSrest = [DYSrest subfreq(:,2,k)]; %#ok<AGROW>
    DYSactive = [DYSactive subfreq(:,4,k)]; %#ok<AGROW>
    
    DYSrestbeta = [DYSrestbeta norm_rest_beta(k)]; %#ok<AGROW>
    DYSactivebeta = [DYSactivebeta norm_active_beta(k)]; %#ok<AGROW>
    
    counterDYS = counterDYS+1;
end

Yrest(:,2) = mean(DYSrest,2);
Erest(:,2) = std(DYSrest,0,2)/(sqrt(counterDYS));
Yactive(:,2) = mean(DYSactive,2);
Eactive(:,2) = std(DYSactive,0,2)/(sqrt(counterDYS));

Ybeta(1,2) = mean(DYSrestbeta);
Ebeta(1,2) = std(DYSrestbeta)/(sqrt(counterDYS));
Ybeta(2,2) = mean(DYSactivebeta);
Ebeta(2,2) = std(DYSactivebeta)/(sqrt(counterDYS));
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
    
    % ETrest and ETactive will contain subfrequency power
    ETrest = [];
    ETactive = [];
    % ETrestbeta and ETactivebeta will contain beta power
    ETrestbeta=[];
    ETactivebeta=[];
    
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
    [norm_rest_beta norm_active_beta subfreq] = ...
        quantPSDshort_loggnorm(rest,active,FREQ_QPSD,order,freq,filename);
        
        %     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
        %now create matrix of percent power in each freq band for each subject
        ETrest = [ETrest subfreq(:,2,k)]; %#ok<AGROW>
        ETactive = [ETactive subfreq(:,4,k)]; %#ok<AGROW>
        
        ETrestbeta = [ETrestbeta norm_rest_beta(k)]; %#ok<AGROW>
        ETactivebeta = [ETactivebeta norm_active_beta(k)]; %#ok<AGROW>
        
        counterET = counterET+1;
    end
    
    Yrest(:,3) = mean(ETrest,2);
    Erest(:,3) = std(ETrest,0,2)/(sqrt(counterET));
    Yactive(:,3) = mean(ETactive,2);
    Eactive(:,3) = std(ETactive,0,2)/(sqrt(counterET));
    
    Ybeta(1,3) = mean(ETrestbeta);
    Ebeta(1,3) = std(ETrestbeta)/(sqrt(counterET));
    Ybeta(2,3) = mean(ETactivebeta);
    Ebeta(2,3) = std(ETactivebeta)/(sqrt(counterET));
end

% SAS 4/29/10: section below commented out
% ------------------------------------------
% %% grab OTHER data
% %Added 4/6/09 for analysis of PD px with dystonic sx
% % if OTHER data needs to be analyzed, change next line to "OTHER=true." Otherwise
% % ET=false;
% OTHER=false;
% if OTHER
% %     pathnameET = uigetdir('','Select directory with OTHER ecogPSD.mat files');
%     OTHERpath = uigetdir('', 'Select directory that contains OTHER _ecogPSD files to be analyzed');
%     OTHERpath = [OTHERpath '\'];
%     cd(OTHERpath);
%     
%     OTHERdir = dir('*_ecogPSD.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script
%     numOTHER = length(OTHERdir);
%     
%     OTHERrest=[];
%     OTHERactive=[];
%     counterOTHER=0;
%     for i=1:numOTHER
%         filename = OTHERdir(i).name;
%         load(filename);
%         if isnan(order(k))
%             continue; %skip files without contact pair over selected area
%         end
%          % for evaluating freq bands other than those already represented in
%     % subfreq matrix, rerun quantPSDmatrix function now, with new freq
%     % bands as defined above
%         [allfreq subfreq] = quantPSDshort_loggnorm(rest,active,FREQ_QPSD,order,freq,filename);
% %     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
%     %now create matrix of percent power in each freq band for each subject
%     OTHERrest = [OTHERrest subfreq(:,2,k)]; %#ok<AGROW>
%         OTHERactive = [OTHERactive subfreq(:,4,k)]; %#ok<AGROW>
%         counterOTHER = counterOTHER+1;
%     end
% 
%     Yrest(:,3) = mean(OTHERrest,2);
% %     Erest(:,3) = std(OTHERrest,0,2);
%     Erest(:,3) = std(OTHERrest,0,2)/sqrt(counterOTHER);
%     Yactive(:,3) = mean(OTHERactive,2);
% %     Eactive(:,3) = std(OTHERactive,0,2);
%     Eactive(:,3) = std(OTHERactive,0,2)/sqrt(counterOTHER);
% end
% ------------------------------------------

%% plot data

% ------------------------FIGURE 1-----------------------
% plot subfrequency power, normalized by gamma power
hf2=figure;
subplot(2,1,1);
handlesr = barEB(Yrest,Erest);
% set(handlesr.ca,'XTick',[1;2;3;4;5],...
%     'XTickLabel',FREQBANDS,...
%     'ylim',YLIM,...
%     'YMinorGrid','on');
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'YMinorGrid','on');
ylim([0 3.5]); % added SAS 4/28/10
title([BRAINAREAS{k} ' at REST']);
ylabel(['Normalized Power' sprintf('\n')...
    '(gamma normalization)']);
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
% set(handlesa.ca,'XTick',[1;2;3;4;5],...
%     'XTickLabel',FREQBANDS,...
%     'ylim',YLIM,...
%     'YMinorGrid','on');
set(handlesa.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'YMinorGrid','on');
ylim([0 3.5]); % added SAS 4/28/10
title([BRAINAREAS{k} ' during ' LIMBAREAS{l} ' MOVEMENT']);
ylabel(['Normalized Power' sprintf('\n')...
    '(gamma normalization)']);
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

% ------------------------FIGURE 1-----------------------
% plot beta  power, normalized by gamma power
hf2=figure;

handlesr = barEB(Ybeta,Ebeta);

set(handlesr.ca,'XTick',[1;2],...
    'XTickLabel',{'REST' 'ACTIVE'},...
    'YMinorGrid','on');
ylim([0 3]); % added SAS 4/28/10
title([BRAINAREAS{k} ' normalized PSD, ' LIMBAREAS{l} ' MOVEMENT']);
ylabel(['Normalized Power' sprintf('\n')...
    'log(beta)-log(gamma)']);
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
% elseif OTHER
%     legend(['PD (n=' num2str(counterPD) ')'],...
%         ['Dys (n=' num2str(counterDYS) ')'],...
%         ['OTHER (n=' num2str(counterOTHER) ')']);
%     hOTHER=handlesr.bars(3);
%     set(hOTHER,'FaceColor','m');
% else
% %     legend(['PD (n=' num2str(counterPD) ')'],...
% %         ['Dys (n=' num2str(counterDYS) ')']);
%     legend(['PD'],...
%         ['Dys']);
end


%% save relevant variables
cd('C:\Users\Starr\Documents\ECOG data\quantPSD data\ANOVA stats');
outputfn=['Ecog_GroupPSDgnorm_' BRAINAREAS{k} '_' LIMBAREAS{l}];
save(outputfn,'PDrest','PDactive','DYSrest','DYSactive','ETrest','ETactive','Yrest','Erest','Yactive','Eactive', 'Ybeta','Ebeta');

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

