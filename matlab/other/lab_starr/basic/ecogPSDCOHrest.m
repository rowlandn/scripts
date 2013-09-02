function ecogPSDCOHrest()
%this program graphs resting state PSD for ecog and lfp channels,
%and ecog-lfp coherence for each channel combination.  It should be run on
%resting activity where there are no movement-related epochs.  It requires
%the output file of apmconv or aomatconv (filename:*_ecog.mat).

% use menu function to input yes or no to whether there is a meaningful lfp
% channel

%note: see mscript runspikeoscil to see robs code on converting timestamps
%to appropriate format (zeros and ones) to run mscohere

%% initialize variables
Fs=1000;                % samples per unit time, in this case
% digitization at 1000 hz rate
% As of 4/21/09, AO data is downsampled from 1.5kHz->1kHz to match sampling
% rates with GL4k data.
WINDOW = 512;           % segment length and Hamming window length for welch's method
NOVERLAP = 256;         % # signal samples that are common to adjacent segments for welch's method
NFFT = 512;             % length of fft
%WINDOW = 4096;           % segment length and Hamming window length for welch's method
%NOVERLAP = 2048;         % # signal samples that are common to adjacent segments for welch's method
%NFFT = 4096;             % length of fft

XLIM_SPEC1 = [0 150];     % range of spectral frequency for plotting (zoomed out)
XLIM_SPEC2 = [0 50];        % range of spectral frequency for plotting (zoomed in)
YLIM_NORM = [0 1.3];    % ylim for plotting normalized PSD
YLIM_LOG = [-3 3];      % ylim for for plotting log PSD
%---frequency ranges for graphing---
FREQ_LO = [8 30];       % beta band
% FREQ_MED = [35 57];     %low gamma band
FREQ_HI = [78 100];     % high gamma band
%---frequency ranges for quantPSD analysis---
% FREQ_QPSD below skips the frequency point at 21.5 and we will include
% that point in low beta for quantPSD analysis (SAS:11/3/09)
% SAS 5/17/10: FREQ_QPSD changed to include 21.5 Hz in lower beta band
% (13-22Hz).
FREQ_QPSD = [4 13;...     % delta alpha band
             13 22;...   % low beta band
             22 31;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band

%% load data
% SS 4/30/09: Instead of using ecog_lfp_raw_data matrix, use ecog structure array
[fn pn]=uigetfile('*_ecog.mat','Select .mat file containing ecog/lfp raw data');
cd(pn);
load([pn fn]);
% remove '_ecog.mat' ending from filename
fn = strrep(fn,'_ecog.mat','');

%% Determine GS4000 or AlphaOmega, set correction factor

system_type = input('Recording system used (1=GS4000, 2=AlphaOmega): ');

if system_type==1 % Guideline 4000 correction factor
    % correction factor based on 1.95Hz frequency resolution
    % (Fs=1000,WINDOW=512) for Guideline 4000 up to 80Hz.
    % Use correction factor below to correct spectral amplitude for frequency
    % points in 1.95-80Hz range
%     crxnfactor = [1e-3,0.0941,0.2936,0.3951,0.4454,0.4926,...
%         0.5888,0.6422,0.6906,0.7454,0.794,0.812,0.8234,...
%         0.8318,0.8389,0.8467,0.8564,0.867,0.8776,0.8876,...
%         0.8964,0.9036,0.9106,0.9174,0.9242,0.9309,0.9374,...
%         0.9437,0.9498,0.9558,0.9614,0.9668,0.9719,0.9767,...
%         0.9811,0.9852,0.9888,0.992,0.9948,0.997,0.9988;];
    % Use correction factor below to only correct spectral amplitude for frequency points above 4Hz
    crxnfactor = [1,1,1,0.3951,0.4454,0.4926,...
        0.5888,0.6422,0.6906,0.7454,0.794,0.812,0.8234,...
        0.8318,0.8389,0.8467,0.8564,0.867,0.8776,0.8876,...
        0.8964,0.9036,0.9106,0.9174,0.9242,0.9309,0.9374,...
        0.9437,0.9498,0.9558,0.9614,0.9668,0.9719,0.9767,...
        0.9811,0.9852,0.9888,0.992,0.9948,0.997,0.9988;];
elseif system_type==2
    % correction factor based on 1.95Hz resolution
    % NOTE: there was no attenuation between 4-100Hz.
    % Use correction factor below to correct spectral amplitude for
    % frequency points in 1.95-80Hz range
%     crxnfactor = [1e-3, 0.791,1,1,1,1,1,1,1,1,1,1,1,1,1,...
%         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    % Use correction factor below to only correct spectral amplitude for
    % frequency points above 4 Hz.
    crxnfactor = ones(1,41);
end
%% re-montage
% ecog12=ecog.contact_pair(1).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal;
% ecog23=ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal;
% ecog34=ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal;
% ecog45=ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(5).raw_ecog_signal;
% ecog56=ecog.contact_pair(5).raw_ecog_signal;
% 
% pt Solis ecog data, input ecog signal and reference contacts were
% accidentally flipped.  Correct for flipped ecog signal amplitude.
%ecog12=ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(1).raw_ecog_signal;
%ecog23=ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal;
%ecog34=ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal;
%ecog45=ecog.contact_pair(5).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal;
%ecog56=-ecog.contact_pair(5).raw_ecog_signal;

% ecogall = [ecog12; ecog23; ecog34; ecog45; ecog56];
% if length(ecog.contact_pair)==6
%     % assume LFP data is in contact pair 6
%     lfp=ecog.contact_pair(6).raw_ecog_signal;
%     ecogall = [ecogall; lfp];
% end

% %-------------------------
% for pt Bujold file 11 who was missing GS4000 Ch.4 (ecog2v6) data
% ecog34=ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal;
% ecog45=ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal;
% ecog56=ecog.contact_pair(4).raw_ecog_signal;
% lfp = ecog.contact_pair(5).raw_ecog_signal;
% ecogall = [ecog34; ecog45; ecog56; lfp];
% %-------------------------

% %-------------------------
% % for pt Nance, 8/21/2009 who had SMA ecog w/ ecog 1 as reference contact
% ecog12=-ecog.contact_pair(2).raw_ecog_signal;
% ecog23=ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal;
% ecog34=ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal;
% ecog45=ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(5).raw_ecog_signal;
% ecog56=ecog.contact_pair(5).raw_ecog_signal-ecog.contact_pair(1).raw_ecog_signal;
% %-------------------------

necog = size(ecogall,1);

%% Define disease state
% Added SAS 3/10/2010
%Will output ecogPSD data into folders specific to disease state
Dx = input('Enter patient diagnosis (1=PD, 2=Dys, 3=ET, 4=Epilepsy, 5=Other): ');
%% Define DBS state during recording 
onoff = input('Was DBS ON during the recording? (y/n) ', 's');

%% PSD using pwelch
%generate power spectral density for each bipolar ecog
%recording and the lfp channel, and coherence spectrum for each ecog channel with lfp
%uses the Welch periodogram with 512 point fft for each

% run pwelch to figure out output length
[psdecog12,freq]=pwelch(ecogall(1,:),WINDOW,NOVERLAP,NFFT,Fs);
psd_length = length(psdecog12);
% initialize psdall
psdall = zeros(psd_length,necog);
for i = 1:necog
    psdall(:,i) = pwelch(ecogall(i,:),WINDOW,NOVERLAP,NFFT,Fs);
    % correct for frequency response5
    for j=1:length(crxnfactor)
        psdall(j,i)=psdall(j,i)/(crxnfactor(j)^2);
    end
end


% normalized PSD
normpsdall = zeros(psd_length,necog); %initialize normpsdall
norm_idx=find(freq>8 & freq<100); % use norm_idx to normalize by max power between 8-100Hz, SAS 11/24/09
for i=1:necog
%     normpsdall(:,i)=psdall(:,i)/max(psdall(:,i)); % normalize each column to its max value
    normpsdall(:,i)=psdall(:,i)/max(psdall(norm_idx(1):norm_idx(end),i)); % normalize each column to its max value
end
% log PSD
logpsdall = log10(psdall);

% %% coherence using mscohere
% 
% % run coherence if LFP exists
% if necog == 6
%     cohall = zeros(psd_length,necog-1);
%     for i = 1:necog-1
%         cohall(:,i) = mscohere(ecogall(i,:),ecogall(6,:), WINDOW,NOVERLAP,NFFT,Fs);
%     end
%     
%     ncoh = size(cohall,2);
%     % transformed coherence
%     tcohall = atanh(sqrt(cohall));
% end

%% quantitative analysis of PSD
% Added SAS 3/10/2010
% now that all the data crunching is done and the results are stored in
% 'psdall' array, run sub-function quantPSD (copied from ecogPSD.m and
% modified for rest data analysis) for quantitative analysis of PSD data

% [allfreq subfreq order] = quantPSD(psdall,FREQ_QPSD,freq,fn,M1_ch);
% 3/12 SAS: Use quantPSDrest.m for quantPSD analysis, not the subfunction above.
[allfreq subfreq order] = quantPSDrest(psdall,FREQ_QPSD,freq,fn,M1_ch);

%% Save and write quantPSD data 
% Added SAS 3/10/2010
% Creates .mat output file of name (filename_ecogPSDrest.mat). 
if Dx==1
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\PD\'];
elseif Dx==2
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\DYS\'];
elseif Dx==3
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\ET\'];
elseif Dx==4
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\Epilepsy\'];
else
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\Other\'];
end
outputname = [outputdir fn,'_ecogPSDrest.mat'];
disp(['Writing ecog PSD data to:  ' outputname]);
% save(outputname,'ecogPSD');
if necog==6
    save(outputname,'psdall', 'tcohall', 'order','freq', 'allfreq', 'subfreq', 'onoff');
else
    save(outputname,'psdall', 'order','freq', 'allfreq', 'subfreq', 'onoff');
end
%% plot data

% first figure
% plot PSD in 4x6 subplot, 6 columns for each ecog/LFP contact, 4 rows for
% different plot views
hf1 = figure;
for i = 1:necog
    
    % 1st row subplot
    subplot(4,necog,i)
    ha = gca;
    hold(ha,'on');
    plot(freq,normpsdall(:,i),'k','LineWidth',2);
    if i==1;
        title([fn sprintf('\n') 'e' num2str(i) num2str(i+1)]); % allows title to have file name
        xlabel('frequency (Hz)');
        ylabel('normalized PSD');
    elseif i == M1_ch
        title(['e' num2str(i) num2str(i+1)],...
            'FontWeight','b');
    elseif i==6
        title('LFP');
    else
        title(['e' num2str(i) num2str(i+1)]);
    end
    set(ha,'YLim',YLIM_NORM);
    set(ha,'XLim',XLIM_SPEC1);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
    
    % 2nd row subplot
    
    subplot(4,necog,necog+i)
    ha = gca;
    hold(ha,'on');
    plot(freq,normpsdall(:,i),'k','LineWidth',2);
    if i==1;
        xlabel('frequency (Hz)');
        ylabel('normalized PSD');
    end
    set(ha,'YLim',YLIM_NORM);
    set(ha,'XLim',XLIM_SPEC2);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
    
    % 3rd row subplot
    subplot(4,necog,2*necog+i)
    ha = gca;
    hold(ha,'on');
    plot(freq,logpsdall(:,i),'k','LineWidth',2);
    if i==1;
        xlabel('frequency (Hz)');
        ylabel('log PSD');
    end
    set(ha,'YLim',YLIM_LOG);
    set(ha,'XLim',XLIM_SPEC1);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
    
    % 4th row subplot
    subplot(4,necog,3*necog+i)
    ha = gca;
    hold(ha,'on');
    plot(freq,logpsdall(:,i),'k','LineWidth',2);
    if i==1;
        xlabel('frequency (Hz)');
        ylabel('log PSD');
    end
    set(ha,'YLim',YLIM_LOG);
    set(ha,'XLim',XLIM_SPEC2);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
end

saveas(hf1,[fn '_ecogPSDrest'],'fig');

% % 2nd figure
% % plot cohere in 4x5 subplot, 5 columns for each contact pair, 4 rows for
% % different plot views
% if necog==6
%     hf2=figure;
%     for i = 1:ncoh
%         % 1st row subplot
%         
%         subplot(4,ncoh,i)
%         ha = gca;
%         hold(ha,'on');
%         plot(freq,cohall(:,i),'k','LineWidth',2);
%         if i==1;
%             title([fn sprintf('\n') 'e' num2str(i) num2str(i+1) '-LFP']); % allows title to have file name
%             xlabel('frequency (Hz)');
%             ylabel('coherence');
%         elseif i==M1_ch
%             title(['e' num2str(i) num2str(i+1) '-LFP'],...
%                 'FontWeight','b');
%         else
%             title(['e' num2str(i) num2str(i+1) '-LFP']);
%         end
%         set(ha,'YLim',[0 1]);
%         set(ha,'XLim',XLIM_SPEC1); % zoomed out
%         FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%         hold(ha,'off');
%         
%         % 2nd row subplot
%         
%         subplot(4,ncoh,ncoh+i)
%         ha = gca;
%         hold(ha,'on');
%         plot(freq,cohall(:,i),'k','LineWidth',2);
%         if i==1
%             xlabel('frequency (Hz)');
%             ylabel('coherence');
%         end
%         set(ha,'YLim',[0 0.5]);
%         set(ha,'XLim',XLIM_SPEC2); % zoomed in
%         FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%         hold(ha,'off');
%         
%         % 3rd row subplot
%         
%         subplot(4,ncoh,2*ncoh+i)
%         ha = gca;
%         hold(ha,'on');
%         plot(freq,tcohall(:,i),'k','LineWidth',2);
%         if i==1
%             xlabel('frequency (Hz)');
%             ylabel('transformed coherence');
%         end
%         set(ha,'YLim',[0 1]);
%         set(ha,'XLim',XLIM_SPEC1); % zoomed in
%         FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%         hold(ha,'off');
%         
%         % 4th row subplot
%         
%         subplot(4,ncoh,3*ncoh+i)
%         ha = gca;
%         hold(ha,'on');
%         plot(freq,tcohall(:,i),'k','LineWidth',2);
%         if i==1
%             xlabel('frequency (Hz)');
%             ylabel('transformed coherence');
%         end
%         set(ha,'YLim',[0 0.5]);
%         set(ha,'XLim',XLIM_SPEC2); % zoomed in
%         FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%         hold(ha,'off');
%         
%     end
%     
%     % save figures
%     saveas(hf2,[fn '_ecogCOHrest'],'fig');
%     
% end
% 
return;

%% sub-functions below
%------------------------------
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

%---------------------------------
function [allfreq subfreq order] = quantPSD(psdall,FREQ_QPSD,freq,filename,M1_ch)
% quantPSD performs quantitative analysis ecogpsdall
% array with the given frequency ranges.  allows user to define ecog
% contact closest to M1
%*UPDATE ALC 6/9/09: no longer allows user to define M1 contact, as this is
%occurs in main function
%
% output: one .mat file with two 3D matrices and # ecog contact that is
% closest to M1.

% select ecog contact closest to M1 and use that to reassign ecog contact
% data array
% ncontact = {'1' '2' '3' '4' '5' '6'};
% contactm1 = menu('Select ecog contact closest to M1',ncontact);
% assign contact numbers for each structure relative to m1 contact
m1 = M1_ch;
% pre = contactm1+1;
% m1s1 = contactm1-2; % updated 2/13:from M1 contact #, subtract 2 instead of 1 to find S1
s1 = m1 - 2; %updated 6/9/09 for clarity

if  m1==1
    s1 = NaN;
    menu(['There is no contact pair over S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
elseif m1==2
    s1 = NaN;
    menu(['There is no contact pair over S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
% elseif contactm1==5
%     pre = NaN;
%     menu(['There is no contact pair over premotor.' sprintf('\n')...
%         'Click OK to continue'],'OK');
elseif m1==6
    m1 = NaN;
%     pre = NaN;
    menu(['There is no contact pair over M1 or premotor.' sprintf('\n')...
        'Click OK to continue'],'OK');
end

% look for LFP channel
num_contact_pair = size(psdall,2);
if num_contact_pair == 6
    lfp = 6; % assume that 6th contact_pair in rest/active structures always contains LFP data
elseif num_contact_pair<6
    lfp = NaN;
end

% order = [pre m1 m1s1 lfp];
order = [m1 s1 lfp];

% initialize and populate the 1x3x4 matrix 'allfreq' (based on ecogPSD.m
% sub-function quantPSD).
% allfreq is a 1x3x4 3-D matrix in which the four 1x3 array corresponds to the
% 4 data recording channels (premotor ecog, M1 ecog, M1-S1 ecog, and stn
% lfp).
% The 1 row contains:
%   1st row     -   resting state data
% The 3 columns contain:
%   1st col     -   frequency at which max power occurs
%   2nd col     -   max power value
%   3rd col     -   total power across all 5 frequency bands

allfreq = zeros(1,3,4);
for i=1:3
    if isnan(order(i))
        allfreq(:,:,i) = NaN;
        continue;
    end
    rest_data = psdall(:,order(i));
% 6/9/09: Need to limit rest  data to those in freq range of interest (typically 4Hz-100Hz)
% 6/19/09: Need to further limit to exclude the frequencies between the bins (ie: exclude 60Hz)    
    rest_data_tmp = rest_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
    freq_tmp = freq(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));    
    [c1 i1] = max(rest_data_tmp);%finds max value of rest data and that value's index position
    allfreq(1,1,i)=freq_tmp(i1);%puts freq corresponding to max value into allfreq matrix
    allfreq(1,2,i)=c1;%puts max value into allfreq matrix
    array1 = [];
    for j = 1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        array1 = [array1; tmp1];  %#ok<AGROW>
    end
    allfreq(1,3,i)=sum(array1);
end

% initialize and populate the 5x2x4 matrix 'subfreq'
% The four 5x2 arrays of subfreq correspond to the 4 data channels
% (pre-motor,M1,M1-S1,STN LFP).  Each of the 5 rows correspond to the 5
% frequency bands defined by variable FREQ_QPSD.
% The 5 columns contain:
%   1st col     -   total power in given frequency band at rest
%   2nd col     -   power in given frequency band, divided by total power
%                   in all 5 freq bands, at rest
subfreq = zeros(5,2,4);
for i=1:3
    if isnan(order(i))
        subfreq(:,:,i)=NaN;
        continue;
    end
    rest_data = psdall(:,order(i));
    for j=1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        subfreq(j,1,i) = sum(tmp1);
        subfreq(j,2,i) = sum(tmp1)/allfreq(1,3,i);
    end
end
% save([filename '_ecogPSD'],'allfreq','subfreq','order');
return;
        


