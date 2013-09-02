%% betadesynch
%graph and define % beta desynchronization across groups
%% define variables for plotting
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   

FREQBANDS = {num2str(FREQ_QPSD(1,:)), num2str(FREQ_QPSD(2,:)) num2str(FREQ_QPSD(3,:)) ...
        num2str(FREQ_QPSD(4,:)) num2str(FREQ_QPSD(5,:))};
% AREAS = {'pre-motor' 'M1' 'S1' 'STN LFP'};
BRAINAREAS = {'M1' 'S1' 'STN LFP'};
LIMBAREAS = {'HAND' 'ELBOW' 'SHOULDER' 'JAW' 'FOOT' '"ARM"' '"NON-ARM"'};
YLIM = [0 70];

%% choose channel
k = menu('Select brain area for analysis',BRAINAREAS);
l = menu('Select limb for analysis',LIMBAREAS);

%% grab PD data
% pathnamePD = uigetdir('','Select directory with PD ecogPSD.mat files');
pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pdpath = [pdpath '\'];
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
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);
    [allfreq subfreq] = quantPSDshort(rest,active,FREQ_QPSD,order,freq,filename);
    %now create matrix of percent power in each freq band for each subject
    PDrest = [PDrest sum(subfreq(2,2,k)+subfreq(3,2,k))];%#ok<AGROW> %2nd column contains % power at rest
    PDactive = [PDactive sum(subfreq(2,4,k)+subfreq(3,4,k))]; %#ok<AGROW> %4th column contains % power with movement
    counterPD = counterPD+1;
end

% PCTdesync = (PDrest - PDactive)./PDrest
ABSdesync = PDactive./PDrest

PDbeta(1,1) = mean(PDrest);
PDbeta(1,2) = mean(PDactive);
PDdesync(1,1) = mean(PDdesync);
%% Plot
h=figure; 