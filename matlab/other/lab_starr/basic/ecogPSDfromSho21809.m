function ecogPSD
% this function calculates PSD using Welch's method to plot the frequency spectrum of files containing digitized
% ecog/LFP data.

% SS 10/2/2008: the montage method assumes that GL4k Ch.3-7 contain ecog data and Ch.8
% contains non-ecog (usually LFP data or noise)

% SS 11/24/2008: updated to perform quantitative analysis of PSD using
% sub-function 'quantPSD.' See 'quantPSD' in the sub-function section below
% for more details.
%% define variables       
Fs=1000;                % samples per unit time, in this case digitization at 1000 hz rate
% Fs = 1.502403855323792e3;   % Alpha Omega ecog/lfp are recorded using different sampling rate. 
% note: in the future, consider downsampling AO data to make it as close to
% GL4k data as possible
OFFSET = 0.5;           % epoch time offset in seconds
EPOCH_LEN = 2;          % epoch length in seconds (1.792 is the lenthg of 6 overlapping segments)
WINDOW = 512;           % segment length and Hamming window length for welch's method
NOVERLAP = 256;         % # signal samples that are common to adjacent segments for welch's method
NFFT = 512;             % length of fft
XLIM_SPEC = [0 150];     % range of spectral frequency for plotting
YLIM_LOG = [-2 2.5];      % ylim for for plotting in log scale
%---frequency ranges for graphing---
FREQ_LO = [8 30];       % beta band
%FREQ_MED = [35 57];     %low gamma band
FREQ_HI = [78 100];     % high gamma band
%---frequency ranges for quantPSD analysis---
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band     

%% import, parse, and montage ecog data
% import filename and pathname of mat file containing ecog data created by APMconv7.m
[fn pn] = uigetfile('*_ecog.mat','Select .mat containing ecog/lfp data');
cd(pn);
% load ecog data
load([pn fn]);
% remove '_ecog.mat' ending from filename
fn = strrep(fn,'_ecog.mat','');

num_contact_pair = length(ecog.contact_pair);
num_epoch = length(ecog.rest_time);

% initialize structures that will contain all rest/active analysis
rest = struct('contact_pair',{});
active = struct('contact_pair',{});
           
for i = 1:num_contact_pair

    % parse each rest/active epoch from each contact pair

    for j = 1:num_epoch
        start_rest = int32((ecog.rest_time(j) + OFFSET) * Fs); % time offset added to epoch times
        end_rest = int32(start_rest + (Fs * EPOCH_LEN) - 1);
        rest(1).contact_pair(i).epoch(j).raw_ecog_signal = ...
            ecog.contact_pair(i).raw_ecog_signal(start_rest:end_rest);
        start_active = int32((ecog.active_time(j)+ OFFSET) * Fs);
        end_active = int32(start_active + Fs*EPOCH_LEN - 1);
        active(1).contact_pair(i).epoch(j).raw_ecog_signal = ...
            ecog.contact_pair(i).raw_ecog_signal(start_active:end_active);
    end

    % montage raw ecog signal to use adjacent contact pairs
    for j = 1: num_epoch
        if i == 1 || i == 6 % the last ecog pair is already 5-6, thus does not require montageing
            break;
        else
            % subtract such that the difference is ecog 1-2, 2-3, 3-4, etc
            % and write over un-montaged raw_ecog_signal field
            rest.contact_pair(i-1).epoch(j).raw_ecog_signal = ...
                (rest.contact_pair(i-1).epoch(j).raw_ecog_signal)...
                - (rest.contact_pair(i).epoch(j).raw_ecog_signal);
            active.contact_pair(i-1).epoch(j).raw_ecog_signal = ...
                (active.contact_pair(i-1).epoch(j).raw_ecog_signal)...
                - (active.contact_pair(i).epoch(j).raw_ecog_signal);
        end
    end

end
%% calculate PSD using pwelch
% perform PSD analysis using pwelch, average all the epochs for
% rest/active

for i = 1:num_contact_pair
    for j = 1:num_epoch  %perform fft for each epoch
        [rest(1).contact_pair(i).epoch(j).PSD ...
            rest(1).contact_pair(i).epoch(j).freq] =...
            pwelch(rest.contact_pair(i).epoch(j).raw_ecog_signal, WINDOW, NOVERLAP, NFFT, Fs);
        [active(1).contact_pair(i).epoch(j).PSD ...
            active(1).contact_pair(i).epoch(j).freq] =...
            pwelch(active.contact_pair(i).epoch(j).raw_ecog_signal,WINDOW,NOVERLAP,NFFT,Fs);
    end

    % calculate mean PSD across all epochs with the subfunction
    % 'mean_signal'
    rest(1).contact_pair(i).mean_PSD = ...
        mean_signal(rest.contact_pair(i),'PSD');
    active(1).contact_pair(i).mean_PSD = ...
        mean_signal(active.contact_pair(i),'PSD');
   
    % normalize mean PSD to peak rest height
    rest(1).contact_pair(i).norm_mean_PSD = ...
        rest.contact_pair(i).mean_PSD/max(rest.contact_pair(i).mean_PSD);
    active(1).contact_pair(i).norm_mean_PSD = ...
        active.contact_pair(i).mean_PSD/max(rest.contact_pair(i).mean_PSD);
   
    % calculate log base 10 of mean PSD
    rest(1).contact_pair(i).log_mean_PSD = ...
        log10(rest.contact_pair(i).mean_PSD);
    active(1).contact_pair(i).log_mean_PSD = ...
        log10(active.contact_pair(i).mean_PSD);
   
    %normalize log mean PSD to peak rest height
    rest(1).contact_pair(i).norm_log_mean_PSD = ...
        rest.contact_pair(i).log_mean_PSD/max(rest.contact_pair(i).log_mean_PSD);
    active(1).contact_pair(i).norm_log_mean_PSD = ...
        active.contact_pair(i).log_mean_PSD/max(rest.contact_pair(i).log_mean_PSD);
   
    % calculate difference in log mean PSD between active and rest
    active.contact_pair(i).difference =...
        active.contact_pair(i).log_mean_PSD - rest.contact_pair(i).log_mean_PSD;
%     active.contact_pair(i).difference =...
%         active.contact_pair(i).norm_log_mean_PSD - rest.contact_pair(i).norm_log_mean_PSD;

    % calculate ratio of active to rest non-normalized mean PSD
%     active.contact_pair(i).ratio =...
%         100*(active.contact_pair(i).mean_PSD./rest.contact_pair(i).mean_PSD);
end

freq = rest(1).contact_pair(1).epoch(1).freq;

%% quantitative analysis of PSD
% now that all the data crunching is done and the results are stored in
% 'rest' and 'active' structure arrays, run sub-function quantPSD for
% quantitative analysis of PSD data
% [allfreq subfreq order] = quantPSD(rest,active,FREQ_QPSD,freq,fn);

%% create figure
% make 5 rows of graphs
% 1st row - plot the averaged PSD, normalized to peak
% 2nd row - same plot as 1st row, expanded 0-50Hz view
% 3rd row - plot the log base 10 of averaged PSD, non-normalized
% 4th row - plot the difference of log averaged PSD
hf1 = figure;
for i = 1:num_contact_pair
   
    % 1st row subplots
    subplot(4,num_contact_pair,i);
    ha = gca;       % create handle for current axes
    hold(ha,'on');
    plot(freq, rest.contact_pair(i).norm_mean_PSD,...
        freq, active.contact_pair(i).norm_mean_PSD,'r','LineWidth',2);
%     plot(freq, rest.contact_pair(i).mean_PSD,...
%         freq, active.contact_pair(i).mean_PSD,'r','LineWidth',2);
    if i == 1
        title([fn sprintf('\n') 'e' num2str(i) '-' num2str(i+1)]);

        ylabel(['normalized average PSD' sprintf('\n')...
            'across all epochs'],...
            'FontSize',8);
        xlabel('frequency (Hz)');
        hl = legend('Rest','Active');
        set(hl,'FontSize',5,'Box','off');
    else
        %         if i == num_contact_pair
        if i == 6 % assume 6th contact_pair contains LFP data
            title('LFP');
        else
            title(['e' num2str(i) '-' num2str(i+1)]);
        end
    end

    set(ha,'YLim',[0 1.3]);
    set(ha,'XLim',XLIM_SPEC);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');

    % 2nd row subplots
    subplot(4,num_contact_pair,num_contact_pair+i);
    ha = gca;       % create handle for current axes
    hold(ha,'on');
    plot(freq, rest.contact_pair(i).norm_mean_PSD,...
        freq, active.contact_pair(i).norm_mean_PSD,'r','LineWidth',2);
%     plot(freq, rest.contact_pair(i).mean_PSD,...
%         freq, active.contact_pair(i).mean_PSD,'r','LineWidth',2);
    if i == 1
        ylabel(['normalized average PSD' sprintf('\n')...
            'across all epochs'],...
            'FontSize',8);
        xlabel('frequency (Hz)');
        hl = legend('Rest','Active');
        set(hl,'FontSize',5,'Box','off');
    end

    set(ha,'YLim',[0 1.3]);
%     set(ha,'YLim',[0 2]); % use for BeatonR27/Jensen
    set(ha,'XLim',[0 50]);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
   
    % 3rd row subplots
    subplot(4,num_contact_pair,2*num_contact_pair+i);
    ha = gca;
    hold(ha,'on');
    plot(freq,rest.contact_pair(i).log_mean_PSD,...
        freq,active.contact_pair(i).log_mean_PSD,...
        'r','LineWidth',2);
%     plot(freq,rest.contact_pair(i).norm_log_mean_PSD,...
%         freq,active.contact_pair(i).norm_log_mean_PSD,...
%         'r','LineWidth',2);
    if i == 1
        ylabel(['log of average' sprintf('\n')...
            'PSD across all epochs'],...
            'FontSize',8);
        xlabel('frequency (Hz)');
        hl = legend('Rest','Active');
        set(hl,'FontSize',5,'Box','off');
    end

    set(ha,'YLim',YLIM_LOG);
    if i==6
        set(ha,'YLim',[-3.5 2.5]); % use for Jensen
    end
    set(ha,'XLim',XLIM_SPEC);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');

    % 4th row subplot
    subplot(4,num_contact_pair,3*num_contact_pair+i);
    ha = gca;
    hold(ha,'on');
    plot(freq, active.contact_pair(i).difference,'r','LineWidth',2)
%    plot(XLIM_SPEC,[100 100],'--k');
   plot(XLIM_SPEC,[0 0],'--k');
    if i == 1
        ylabel(['difference of' sprintf('\n')...
            'log average PSD' sprintf('\n')...
            'between active and rest'],...
            'FontSize',8);
        xlabel('frequency (Hz)');
        hl = legend('Active');
        set(hl,'FontSize',5,'Box','off');
    end
    set(ha,'XLim',XLIM_SPEC);
%     set(ha,'YLim',[0 800]); % keep the same y scale for all plots to facilitate comparison across contacts
    set(ha,'YLim',[-2 1]); % keep the same y scale for all plots to
    % facilitate comparison across contacts
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma bands
    hold(ha,'off');

end


% save figure
saveas(hf1,[fn '_ecogPSD'],'fig');

%% create optional figures
% %this plot visualizes data in allfreq
% hf2 = figure;
% areas = {'pre-motor' 'M1' 'M1-S1' 'LFP'};
% freqbands = {'4-12' '13-21' '22-30' '31-55' '76-100'};
% for i=1:4
%     subplot(4,1,i);
%     if isnan(order(i));
%         text(0.5,0.5,...
%             [areas{i} ' not available'],...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle');
%         continue;
%     end
%     ha = gca;
%     hold(ha,'on');
%     plot(freq, rest.contact_pair(order(i)).mean_PSD,...
%         freq, active.contact_pair(order(i)).mean_PSD,'r','LineWidth',2);
%     plot(allfreq(1,1,i),allfreq(1,2,i),'o');
%     text(allfreq(1,1,i),allfreq(1,2,i),['Max Freq: ' num2str(allfreq(1,1,i)) sprintf('\n')...
%         'Power Value: ' num2str(allfreq(1,2,i))]);
%     plot(allfreq(2,1,i),allfreq(2,2,i),'or');
%     text(allfreq(2,1,i),allfreq(2,2,i),['Max Freq: ' num2str(allfreq(2,1,i)) sprintf('\n')...
%         'Power Value: ' num2str(allfreq(2,2,i))]);
%
%     set(ha,'XLim',[0 50]);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     title(areas{i});
%     hold(ha,'off');
% end
%
% %this plot visualizes data in subfreq
% hf3 = figure;
% x = [1 2 3 4 5];
% for i=1:4
%     subplot(4,1,i);
%     if isnan(order(i));
%         text(0.5,0.5,...
%             [areas{i} ' not available'],...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle');
%         continue;
%     end
%     ha=gca;
%     hold(ha,'on');
%     Y(:,1)=subfreq(:,2,i);
%     Y(:,2)=subfreq(:,4,i);
%     bar(x,Y);
%     set(ha,'XTick',[1;2;3;4;5]);
%     set(ha,'XTickLabel',freqbands);
%     title(areas{i});
%     hold(ha,'off');
% end


%% sub-functions below
function out = mean_signal(contact_pair, fieldname)
num_epoch = length(contact_pair.epoch);
A = [];
for j = 1:num_epoch
    tmp = eval(['contact_pair.epoch(' num2str(j) ').' fieldname])';
    A = [A;tmp]; %#ok<AGROW> % concatenate
end
out = mean(A);
return;

function FFT_band_fill(ha,FREQ_HI,FREQ_LO)
% FFT_band_fill fills bands of specified frequncies on axes specified by
% the handle ha

ylm = get(ha,'YLim');
x_lo = [FREQ_LO(1) FREQ_LO(2) FREQ_LO(2) FREQ_LO(1)];
x_hi = [FREQ_HI(1) FREQ_HI(2) FREQ_HI(2) FREQ_HI(1)];
y = [ylm(1) ylm(1) ylm(2) ylm(2)];
fill(x_lo,y,'g',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
fill(x_hi,y,'y',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
return;

function [allfreq subfreq order] = quantPSD(rest,active,FREQ_QPSD,freq,filename)
% quantPSD performs quantitative analysis on rest/active structure
% arrays with the given frequency ranges.  allows user to define ecog
% contact closest to M1
%
% output: one .mat file with two 3D matrices and # ecog contact that is
% closest to M1.

% select ecog contact closest to M1 and use that to reassign ecog contact
% data array
ncontact = {'1' '2' '3' '4' '5' '6'};
contactm1 = menu('Select ecog contact closest to M1',ncontact);
% assign contact numbers for each structure relative to m1 contact
m1 = contactm1;
pre = contactm1+1;
m1s1 = contactm1-2; % updated 2/13:from M1 contact #, subtract 2 instead of 1 to find S1

if contactm1==1
    m1s1 = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
elseif contactm1==5
    pre = NaN;
    menu(['There is no contact pair over premotor.' sprintf('\n')...
        'Click OK to continue'],'OK');
elseif contactm1==6
    m1 = NaN;
    pre = NaN;
    menu(['There is no contact pair over M1 or premotor.' sprintf('\n')...
        'Click OK to continue'],'OK');
end

% look for LFP channel
num_contact_pair = length(rest.contact_pair);
if num_contact_pair == 6
    lfp = 6; % assume that 6th contact_pair in rest/active structures always contains LFP data
elseif num_contact_pair<6
    lfp = NaN;
end

order = [pre m1 m1s1 lfp];

% initialize and populate the 2x3x4 matrix 'allfreq'
% allfreq is a 2x3x4 3-D matrix in which the four 2x3 arrays correspond to the
% 4 data recording channels (premotor ecog, M1 ecog, M1-S1 ecog, and stn
% lfp).
% The 2 rows contain:
%   1st row     -   resting state data
%   2nd row     -   active movement data
% The 3 columns contain:
%   1st col     -   frequency at which max power occurs
%   2nd col     -   max power value
%   3rd col     -   total power across all 5 frequency bands

allfreq = zeros(2,3,4);
for i=1:4
    if isnan(order(i))
        allfreq(:,:,i) = NaN;
        continue;
    end
    rest_data = rest.contact_pair(order(i)).mean_PSD;
    active_data = active.contact_pair(order(i)).mean_PSD;
    [c1 i1] = max(rest_data);
    allfreq(1,1,i)=freq(i1);
    allfreq(1,2,i)=c1;
    [c2 i2] = max(active_data);
    allfreq(2,1,i)=freq(i2);
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
subfreq = zeros(5,5,4);
for i=1:4
    if isnan(order(i))
        subfreq(:,:,i)=NaN;
        continue;
    end
    rest_data = rest.contact_pair(order(i)).mean_PSD;
    active_data = active.contact_pair(order(i)).mean_PSD;
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
save([filename '_ecogPSD'],'allfreq','subfreq','order');
return;

   





-----Inline Attachment Follows-----

function PlotQuantPSD
% plot results of quantPSD output from ecogPSD analysis

%%define variables
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   
FREQBANDS = {'4-12' '13-21' '22-30' '31-55' '76-100'};
AREAS = {'pre-motor' 'M1' 'M1-S1' 'STN LFP'};
YLIM = [0 70];
%% choose channel
k = menu('Select area for analysis',AREAS);

%% grab dystonia data
% pathnameDYT = uigetdir('','Select directory with dystonia ecogPSD.mat files');
pathnameDYT = 'C:\Documents and Settings\Sho\My Documents\Lab Documents\Ecog-LFP-SSEP\DMRF grant\dystonia bi elbow';
cd(pathnameDYT);
FileListDYT = dir('*_ecogPSD.mat');
nfilesDYT = length(FileListDYT);
DYTrest=[];
DYTactive=[];
counterDYT=0;
for i=1:nfilesDYT
    filename = FileListDYT(i).name;
    load(filename);
    if isnan(order(k))
        continue; %skip files without contact pair over selected area
    end
    DYTrest = [DYTrest subfreq(:,2,k)]; %#ok<AGROW>
    DYTactive = [DYTactive subfreq(:,4,k)]; %#ok<AGROW>
    counterDYT = counterDYT+1;
end

Yrest(:,1) = mean(DYTrest,2);
Erest(:,1)= std(DYTrest,0,2);
Yactive(:,1) = mean(DYTactive,2);
Eactive(:,1) = std(DYTactive,0,2);

%% grab PD data
% pathnamePD = uigetdir('','Select directory with PD ecogPSD.mat files');
pathnamePD='C:\Documents and Settings\Sho\My Documents\Lab Documents\Ecog-LFP-SSEP\DMRF grant\PD elbow';
cd(pathnamePD);
FileListPD = dir('*_ecogPSD.mat');
nfilesPD = length(FileListPD);
PDrest = [];
PDactive=[];
counterPD = 0;
for i=1:nfilesPD
    filename = FileListPD(i).name;
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
ET=true;
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
    legend(['Dyst (n=' num2str(counterDYT) ')'],...
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





Go to Previous message | Go to Next message | Back to Messages
ASCII (ASCII)Greek (ISO-8859-7)Greek (Windows-1253)Latin-10 (ISO-8859-16)Latin-3 (ISO-8859-3)Latin-6 (ISO-8859-10)Latin-7 (ISO-8859-13)Latin-8 (ISO-8859-14)Latin-9 (ISO-8859-15)W. European (850)W. European (CP858)W. European (HPROMAN8)W. European (MACROMAN8)W. European (Windows-1252)Armenia (ARMSCII-8)Baltic Rim (ISO-8859-4)Baltic Rim (WINDOWS-1257)Cyrillic (866)Cyrillic (ISO-8859-5)Cyrillic (KOI8-R)Cyrillic (KOI8-RU)Cyrillic (KOI8-T)Cyrillic (KOI8-U)Cyrillic (WINDOWS-1251)Latin-2 (852)Latin-2 (ISO-8859-2)Latin-2 (WINDOWS-1250)Turkish (ISO-8859-9)Turkish (WINDOWS-1254)Arabic (ISO-8859-6, ASMO-708)Arabic (WINDOWS-1256)Hebrew (856)Hebrew (862)Hebrew (WINDOWS-1255)Chinese Simplified (GB-2312-80)Chinese Simplified (GB18030)Chinese Simplified (HZ-GB-2312)Chinese Simplified (ISO-2022-CN)Chinese Simplified (WINDOWS-936)Chinese Trad.-Hong Kong (BIG5-HKSCS)Chinese Traditional (BIG5)Chinese Traditional (EUC-TW)Japanese (SHIFT_JIS)Japanese (EUC-JP)Japanese (ISO-2022-JP)Korean (ISO-2022-KR)Korean (EUC-KR)Thai (TIS-620-2533)Thai (WINDOWS-874)Vietnamese (TCVN-5712)Vietnamese (VISCII)Vietnamese (WINDOWS-1258)Unicode (UTF-7)Unicode (UTF-8)Unicode (UTF-16)Unicode (UTF-32)
| Full Headers

Search Mail


Copyright � 1994-2009 Yahoo! Inc. All rights reserved. Terms of Service - Copyright/IP Policy - Guidelines
NOTICE: We collect personal information on this site.
To learn more about ho