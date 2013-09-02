function Ecog_GroupPSDabs
% Definition: same as Ecog_GroupPSD, except bars represent the actual power
% in each frequency band, rather than the percent power relative to the
% total power across the five frequency bands - ALC 3/9/10

% INPUT = group _ecogPSD files from quantPSD data folder
% specifically, needs mean PSD data

% CALLS function 'quantPSDshort' - this is nearly identical to the
% quantPSD subfunction called by ecogPSD, except that it does not assign
% contact numbers or create the variable 'order' - as it already exists.
% The purpose of calling 'quantPSDshort'is to recalcuate the allfreq and
% subfreq variables in the event that the user wishes to change the
% FREQ_QPSD variable for the purpose of analyzing group data.
% 'quantPSDshort' also does not write data to an excel file.

% OUTPUT: mat file with variables for ANOVA analysis, creates plot

% 3/9/10: this is identical to Ecog_GroupPSD except instead of graphing the
% % power of each freq band (subfreq(:,2) for rest and subfreq(:,4)) for
% active state, it is looking for "absolute" power (subfreq(:,1) and
% subfreq(:,3) respectively)**ALC
%% define variables for plotting
% SAS 4/29/10: changed frequency bandwidth so that they are contiguous (ie.
% 4-13,13-22,22-31,31-55,76-100).  Note:21.48Hz point is included in low
% beta, which gives 5 total data points in low gamma at frequency
% resolution of 1.95 Hz.
FREQ_QPSD = [4 13;...     % delta alpha band
             13 20;...   % low beta band
             20 30;...   % high beta band
             30 55;...   % low gamma band
             76 100];    % high gamma band   
%          FREQ_QPSD = [136 160;...     % delta alpha band
%              196 220;...   % low beta band
%              256 280;...   % high beta band
%              316 340;...   % low gamma band
%              436 460];    % high gamma band  
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
% pdpath = 'C:\Users\Elena\Documents\ECOG data\quantPSD data\PD\PSD corrected for Andrea''s study';
pdpath = 'C:\Users\Elena\Desktop\Random stuff for Andrea''s paper\Original data\PD elbow';
% pdpath = 'C:\Users\Elena\Desktop\Gamma 1 sec epoch analysis\PD';
% pdpath = [pdpath '\' LIMBAREAS{l}];
cd(pdpath);

PDdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numPD = length(PDdir);

PDrest = [];
PDactive=[];
PDdiff=[];
PDbetarest = [];
PDbetaactive = [];
PDgammarest = [];
PDgammaactive = [];
PDdiff1 = [];
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
order(1) = order(1);
    if make_logplot
        [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename);
    else
        [allfreq subfreq] = quantPSDshort(rest,active,FREQ_QPSD,order,freq,filename);
    end
    
    %now create matrix of percent power in each freq band for each subject
    %to calculate "absolute" power in each band (rather than "% power") use
    %subfreq(:,1,k) and subfreq(:,3,k) for rest and active respectively
    PDrest = [PDrest subfreq(:,1,k)];%#ok<AGROW> %1st column contains absolute power at rest
    PDactive = [PDactive subfreq(:,3,k)]; %#ok<AGROW> %3rd column contains absolute power with movement
    PDdiff = [PDdiff subfreq(:,5,k)];
    PDdiff1 = [PDdiff1 subfreq(:,1,k)-subfreq(:,3,k)];
    PDbetarest = [PDbetarest (subfreq(2,1,k) + subfreq(3,1,k))];
    PDbetaactive = [PDbetaactive (subfreq(2,3,k) + subfreq(3,3,k))];
    PDgammarest = [PDgammarest subfreq(5,1,k)];
    PDgammaactive = [PDgammaactive subfreq(5,3,k)];
    counterPD = counterPD+1;
    
end


Yrest(:,1) = mean(PDrest,2);
Erest(:,1) = std(PDrest,0,2)/(sqrt(counterPD));
RestBetaPower(:,1) = PDbetarest';
RestBetaMean = mean(PDbetarest);
RestBetaMeanE = std(PDbetarest)/(sqrt(counterPD));

Yactive(:,1) = mean(PDactive,2);
Eactive(:,1) = std(PDactive,0,2)/(sqrt(counterPD));
ActiveBetaPower(:,1) = PDbetaactive';
ActiveBetaMean = mean(PDbetaactive);
ActiveBetaMeanE = std(PDbetaactive)/(sqrt(counterPD));


%Ydiff is a 5x3 matrix.  5 rows for 5 frequency bands.  3 columns for PD,
%DYST, and ET
Ydiff(:,1)= mean(PDdiff,2);
Ediff(:,1)= std(PDdiff,0,2)/sqrt(counterPD);
% Ybeta is a 2x3 matrix.  2 rows for rest and active.  3 columns for PD, DYS,
% and ET
Ybeta(1,1) = mean(PDbetarest);
Ebeta(1,1) = std(PDbetarest)/(sqrt(counterPD));
Ybeta(2,1) = mean(PDbetaactive);
Ebeta(2,1) = std(PDbetaactive)/(sqrt(counterPD));

Ygamma(1,1) = mean(PDgammarest);
Egamma(1,1) = std(PDgammarest)/(sqrt(counterPD));
Ygamma(2,1) = mean(PDgammaactive);
Egamma(2,1) = std(PDgammaactive)/(sqrt(counterPD));


%% grab dystonia data
% pathnameDYT = uigetdir('','Select directory with dystonia ecogPSD.mat files');
% dyspath = uigetdir('', 'Select directory that contains Dys _ecgPSD files to be analyzed');
% dyspath = 'C:\Users\Elena\Documents\ECOG data\quantPSD data\DYS\PSD corrected for Andrea''s study';
dyspath = 'C:\Users\Elena\Desktop\Random stuff for Andrea''s paper\Original data\DYS elbow';
% dyspath = 'C:\Users\Elena\Desktop\Gamma 1 sec epoch analysis\DYS';
% dyspath = [dyspath '\' LIMBAREAS{l}];
cd(dyspath);

DYSdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
% FileListDYS = dir('*_ecogPSD.mat');

% nfilesDYT = length(FileListDYT);
numDYS = length(DYSdir);

DYSrest=[];
DYSactive=[];
DYSdiff=[];
DYSbetarest = [];
DYSbetaactive = [];
DYSgammaactive = [];
DYSgammarest = [];


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
order(1) = order(1);
    if make_logplot
        [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename);
    else
        [allfreq subfreq] = quantPSDshort(rest,active,FREQ_QPSD,order,freq,filename);
    end
    %now create matrix of percent power in each freq band for each subject
    %to calculate "absolute" power in each band (rather than "% power") use
    %subfreq(:,1,k) and subfreq(:,3,k) for rest and active respectively
    DYSrest = [DYSrest subfreq(:,1,k)]; %#ok<AGROW>
    DYSactive = [DYSactive subfreq(:,3,k)]; %#ok<AGROW>
    DYSdiff = [DYSdiff subfreq(:,5,k)]; %#ok<AGROW>
    DYSbetarest = [DYSbetarest (subfreq(2,1,k) + subfreq(3,1,k))];
    DYSbetaactive = [DYSbetaactive (subfreq(2,3,k) + subfreq(3,3,k))];
    DYSgammarest = [DYSgammarest subfreq(5,1,k)];
    DYSgammaactive = [DYSgammaactive subfreq(5,3,k)];
    counterDYS = counterDYS+1;
    
end

Yrest(:,2) = mean(DYSrest,2);
Erest(:,2) = std(DYSrest,0,2)/(sqrt(counterDYS));
%RestBetaPower(:,2) = DYSbetarest';
RestBetaMean = [RestBetaMean mean(DYSbetarest)];
RestBetaMeanE = [RestBetaMeanE std(DYSbetarest)/(sqrt(counterDYS))];

Yactive(:,2) = mean(DYSactive,2);
Eactive(:,2) = std(DYSactive,0,2)/(sqrt(counterDYS));
%ActiveBetaPower(:,2) = DYSbetaactive';
ActiveBetaMean = [ActiveBetaMean mean(DYSbetaactive)];
ActiveBetaMeanE = [ActiveBetaMeanE std(DYSbetaactive)/(sqrt(counterDYS))];

Ydiff(:,2)= mean(DYSdiff,2);
Ediff(:,2)= std(DYSdiff,0,2)/sqrt(counterDYS);

Ybeta(1,2) = mean(DYSbetarest);
Ebeta(1,2) = std(DYSbetarest)/(sqrt(counterDYS));
Ybeta(2,2) = mean(DYSbetaactive);
Ebeta(2,2) = std(DYSbetaactive)/(sqrt(counterDYS));

Ygamma(1,2) = mean(DYSgammarest);
Egamma(1,2) = std(DYSgammarest)/(sqrt(counterDYS));
Ygamma(2,2) = mean(DYSgammaactive);
Egamma(2,2) = std(DYSgammaactive)/(sqrt(counterDYS));

%% grab ET data
% if ET data needs to be analyzed, change next line to "ET=true." Otherwise
% ET=false;
ET=true;
if ET
    %     pathnameET = uigetdir('','Select directory with ET ecogPSD.mat files');
    %     etpath = uigetdir('', 'Select directory that contains ET _ecogPSD files to be analyzed');
%     etpath = 'C:\Users\Elena\Documents\ECOG data\quantPSD data\ET\corrected PSD for Andrea''s study';
    etpath = 'C:\Users\Elena\Desktop\Random stuff for Andrea''s paper\Original data\ET elbow';
% etpath = 'C:\Users\Elena\Desktop\Gamma 1 sec epoch analysis\ET';
%     etpath = [etpath '\' LIMBAREAS{l}];
    cd(etpath);
    
    ETdir = dir('*_ecogPSD.mat'); % selects '*_transcoh.mat files' that were outputs of the coherenece analysis script
    numET = length(ETdir);
    
    ETrest=[];
    ETactive=[];
    ETdiff=[];
    ETbetarest = [];
    ETbetaactive = [];
    ETgammarest = [];
    ETgammaactive = [];
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
        order(1) = order(1);
        if make_logplot
            [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename);
        else
            [allfreq subfreq] = quantPSDshort(rest,active,FREQ_QPSD,order,freq,filename);
        end
        %     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
        %now create matrix of percent power in each freq band for each subject
        %to calculate "absolute" power in each band (rather than "% power") use
        %subfreq(:,1,k) and subfreq(:,3,k) for rest and active respectively
        ETrest = [ETrest subfreq(:,1,k)]; %#ok<AGROW>
        ETactive = [ETactive subfreq(:,3,k)]; %#ok<AGROW>
        ETdiff = [ETdiff subfreq(:,5,k)]; %#ok<AGROW>
        ETbetarest = [ETbetarest (subfreq(2,1,k) + subfreq(3,1,k))];
        ETbetaactive = [ETbetaactive (subfreq(2,3,k) + subfreq(3,3,k))];
        ETgammarest = [ETgammarest subfreq(5,1,k)];
        ETgammaactive = [ETgammaactive subfreq(5,3,k)];
        counterET = counterET+1;
    end
    
    Yrest(:,3) = mean(ETrest,2);
    Erest(:,3) = std(ETrest,0,2)/(sqrt(counterET));
%    RestBetaPower(:,3) = ETbetarest';
    RestBetaMean = [RestBetaMean mean(ETbetarest)];
    RestBetaMeanE = [RestBetaMeanE std(ETbetarest)/(sqrt(counterET))];
    
    Yactive(:,3) = mean(ETactive,2);
    Eactive(:,3) = std(ETactive,0,2)/(sqrt(counterET));
    %ActiveBetaPower(:,1) = ETbetaactive';
    ActiveBetaMean = [ActiveBetaMean mean(ETbetaactive)];
    ActiveBetaMeanE = [ActiveBetaMeanE std(ETbetaactive)/(sqrt(counterET))];
    
    Ydiff(:,3) = mean(ETdiff,2);
    Ediff(:,3) = std(ETdiff,0,2)/(sqrt(counterET));
    
    Ybeta(1,3) = mean(ETbetarest);
    Ebeta(1,3) = std(ETbetarest)/(sqrt(counterET));
    Ybeta(2,3) = mean(ETbetaactive);
    Ebeta(2,3) = std(ETbetaactive)/(sqrt(counterET));
    
    Ygamma(1,3) = mean(ETgammarest);
    Egamma(1,3) = std(ETgammarest)/(sqrt(counterET));
    Ygamma(2,3) = mean(ETgammaactive);
    Egamma(2,3) = std(ETgammaactive)/(sqrt(counterET));
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
        if make_logplot
            [allfreq subfreq] = quantPSDshort_log(rest,active,FREQ_QPSD,order,freq,filename);
        else
            [allfreq subfreq] = quantPSDshort(rest,active,FREQ_QPSD,order,freq,filename);
        end
        %     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
        %now create matrix of percent power in each freq band for each subject
        OTHERrest = [OTHERrest subfreq(:,1,k)]; %#ok<AGROW>
        OTHERactive = [OTHERactive subfreq(:,3,k)]; %#ok<AGROW>
        counterOTHER = counterOTHER+1;
    end
    
    Yrest(:,3) = mean(OTHERrest,2);
    Erest(:,3) = std(OTHERrest,0,2)/sqrt(counterOTHER);
    
    Yactive(:,3) = mean(OTHERactive,2);
    Eactive(:,3) = std(OTHERactive,0,2)/sqrt(counterOTHER);
    
    Ydiff(:,3) = mean(ETdiff,2);
    Ediff(:,3) = std(ETdiff,0,2)/(sqrt(counterET));
end

%% plot data
% ------------------FIGURE 1--------------------------
% plot rest/active abs log power
hf1=figure;
%temporary code for graphing only beta band
% Yrest= Yrest(2,:);
% Erest = Erest(2,:); 
% %------
subplot(2,1,1);
handlesr = barEB(Yrest,Erest);
% set(handlesr.ca,'XTick',[1;2;3;4;5],...
%     'XTickLabel',FREQBANDS,...
%     'ylim',YLIM,...
%     'YMinorGrid','on');
set(gca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'YMinorGrid','on');
if make_logplot
    ylim([-2 2.5]); % added SAS 4/28/10
end
title([BRAINAREAS{k} ' at REST']);
ylabel('Power');
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
set(gca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'YMinorGrid','on');
if make_logplot
    ylim([-2 2.5]); % added SAS 4/28/10
end
title([BRAINAREAS{k} ' during ' LIMBAREAS{l} ' MOVEMENT']);
ylabel('Power');
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

%-----------------FIGURE 2--------------------
% plot movement related power change in log scale
% plot rest/active abs log power
hf2=figure;
handlesr = barEB(Ydiff,Ediff);
set(gca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'YMinorGrid','on');

if make_logplot
    ylim([-1.5 0.5]); % added SAS 4/28/10
end

title([BRAINAREAS{k} ', Movement-related power change' sprintf('\n')...
    'during ' LIMBAREAS{l} ' MOVEMENT']);
ylabel(['Power change' sprintf('\n')...
    '(ACTIVE - REST)']);
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
%----------------------Figure 3------------------------------
%hf3=figure;
%temporary code for graphing only beta band
% Yrest= Yrest(2,:);
% Erest = Erest(2,:); 
% %------
%handlesr = barEB(Ybeta,Ebeta);
% set(handlesr.ca,'XTick',[1;2;3;4;5],...
%     'XTickLabel',FREQBANDS,...
%     'ylim',YLIM,...
%     'YMinorGrid','on');
%set(gca,'XTick',[1;2],...
   %'XTickLabel',{'REST','ACTIVE'},...
   %'YMinorGrid','on');
%if make_logplot
   % ylim([-2 6]); % added SAS 4/28/10
%end
%title([BRAINAREAS{k} ' total beta power']);
%ylabel('Power');

%hPD=handlesr.bars(1);
%hDYS=handlesr.bars(2);
%set(hPD,'FaceColor','b');
%set(hDYS,'FaceColor','g');
%if ET
    %legend(['PD (n=' num2str(counterPD) ')'],...
        %['Dys (n=' num2str(counterDYS) ')'],...
        %['ET (n=' num2str(counterET) ')']);
   % hET=handlesr.bars(3);
    %set(hET,'FaceColor','m');
%elseif OTHER
   % legend(['PD (n=' num2str(counterPD) ')'],...
        %['Dys (n=' num2str(counterDYS) ')'],...
        %['OTHER (n=' num2str(counterOTHER) ')']);
   %hOTHER=handlesr.bars(3);
   % set(hOTHER,'FaceColor','m');
%else
%     legend(['PD (n=' num2str(counterPD) ')'],...
%         ['Dys (n=' num2str(counterDYS) ')']);
    %legend(['PD'],...
        %['Dys']);
%end

%% save relevant variables
cd('C:\Users\Elena\Documents\ECOG data\quantPSD data\ANOVA stats');
outputfn=['Ecog_GroupPSDabs_' BRAINAREAS{k} '_' LIMBAREAS{l}];
save(outputfn,'PDrest','PDactive','DYSrest','DYSactive','ETrest','ETactive','Yrest','Erest','Yactive','Eactive', 'Ydiff','Ediff');

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

