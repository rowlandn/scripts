function ecogPSD_rest_highgrid(name)
%this program graphs resting state PSD for ecog and lfp channels,
%and ecog-lfp coherence for each channel combination.  It should be run on
%resting activity where there are no movement-related epochs.  It requires
%the output file of apmconv or aomatconv (filename:*_ecog.mat).

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

    % ylim for for plotting log PSD
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
% [fn pn]=uigetfile('*_ecog.mat','Select .mat file containing ecog/lfp raw data');
% cd(pn);
% load([pn fn]);
% remove '_ecog.mat' ending from filename
 load(name);
fn = strrep(name,'_ecog.mat','');

necog = length(ecog.contact_pair);
 %% Define disease state
% % Added SAS 3/10/2010
% %Will output ecogPSD data into folders specific to disease state
% Dx = input('Enter patient diagnosis (1=PD, 2=Dys, 3=ET, 4=Epilepsy, 5=Other): ');
% %% Define DBS state during recording 
% onoff = input('Was DBS ON during the recording? (y/n) ', 's');

%% PSD using pwelch
%generate power spectral density for each bipolar ecog
%recording and the lfp channel, and coherence spectrum for each ecog channel with lfp
%uses the Welch periodogram with 512 point fft for each

% run pwelch to figure out output length
[psdecog12,freq]=pwelch(ecog.contact_pair(1).remontaged_ecog_signal,WINDOW,NOVERLAP,NFFT,Fs);
psd_length = length(psdecog12);
% initialize psdall
psdall = zeros(psd_length,necog);
for i = 1:necog
    psdall(:,i) = pwelch(ecog.contact_pair(i).remontaged_ecog_signal,WINDOW,NOVERLAP,NFFT,Fs);
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

%% plot data

XLIM_SPEC1 = [0 150];     % range of spectral frequency for plotting (zoomed out)
XLIM_SPEC2 = [0 50];% range of spectral frequency for plotting (zoomed in)
x= max(max(psdall));
YLIM_RAW = [0 x];    % ylim for plotting raw PSD
x = max(max(normpsdall));
YLIM_NORM = [0 x];    % ylim for plotting normalized PSD
x = max(max(logpsdall));
YLIM_LOG = [-x x];  
if necog ==5
   length_row = 5;
elseif necog ==6
   length_row = 6;
elseif necog <=32
    length_row = 14;
    M1_ch=M1_ch1;
elseif necog >= 64;
    length_row = 16;
    M1_ch=M1_ch1;
end

% % figure plot raw PSD 
% row=1;
% for ii = 1:length_row:necog
%     hf1 = figure;
%     for i = 1:length_row
%         chan=ii+i-1;
% %         M1_ch=M1_ch1;
%         % 1st row subplot
%         subplot(4,4,i)
%         ha = gca;
%         hold(ha,'on');
%         plot(freq,psdall(:,chan),'k','LineWidth',2);
%         if i==1;
%             title([fn sprintf('\n') 'e' num2str(chan) ]); % allows title to have file name
%             xlabel('frequency (Hz)');
%             ylabel('raw PSD');
%         elseif i == M1_ch
%             title(['e' num2str(chan) ],...
%                 'FontWeight','b');
%         else
%             title(['e' num2str(chan) ]);
%         end
%         set(ha,'YLim',YLIM_RAW);
%         set(ha,'XLim',XLIM_SPEC1);
%         FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%         hold(ha,'off');
%     end
%     row = row;
%     saveas(hf1,[fn '_rawpsd' num2str(ii)],'fig');
%     row = row+1;
% end
% 
% % figure plot normalized PSD 
% row=1;
% for ii = 1:length_row:necog
%     hf1 = figure;
%     for i = 1:length_row
%         chan=ii+i-1;
% %         M1_ch=M1_ch1;
%         % 1st row subplot
%         subplot(4,4,i)
%         ha = gca;
%         hold(ha,'on');
%         plot(freq,normpsdall(:,chan),'k','LineWidth',2);
%         if i==1;
%             title([fn sprintf('\n') 'e' num2str(chan) ]); % allows title to have file name
%             xlabel('frequency (Hz)');
%             ylabel('norm PSD');
%         elseif i == M1_ch
%             title(['e' num2str(chan) ],...
%                 'FontWeight','b');
%         else
%             title(['e' num2str(chan) ]);
%         end
%         set(ha,'YLim',YLIM_NORM);
%         set(ha,'XLim',XLIM_SPEC1);
%         FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%         hold(ha,'off');
%     end
%     row = row;
%     saveas(hf1,[fn '_normpsd' num2str(ii)],'fig');
%     row = row+1;
% end

% figure plot log PSD 
row=1;
for ii = 1:length_row:necog
    hf1 = figure;
    for i = 1:length_row
        chan=ii+i-1;
%         M1_ch=M1_ch1;
        % 1st row subplot
        subplot(4,4,i)
        ha = gca;
        hold(ha,'on');
        plot(freq,logpsdall(:,chan),'k','LineWidth',2);
        if i==1;
            title([fn sprintf('\n') 'e' num2str(chan) ]); % allows title to have file name
            xlabel('frequency (Hz)');
            ylabel('log PSD');
        elseif i == M1_ch
            title(['e' num2str(chan) ],...
                'FontWeight','b');
        else
            title(['e' num2str(chan) ]);
        end
        set(ha,'YLim',YLIM_LOG);
        set(ha,'XLim',XLIM_SPEC1);
        FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
        hold(ha,'off');
    end
    row = row;
    saveas(hf1,[fn '_logpsd' num2str(ii)],'fig');
    row = row+1;
end

close all
% % third figure plot PSD in 4x6 subplot, 6 columns for ecog/LFP contact around M1, 4 rows for
% % different plot views
% hf1 = figure;
% M1=M1_ch1;
% x=[M1-2 M1-1 M1 M1+1 M1+2];
% for i = x
%     j=i-x(1)+1;
%     % 1st row subplot
%     subplot(4,5,j)
%     ha = gca;
%     hold(ha,'on');
%     plot(freq,normpsdall(:,i),'k','LineWidth',2);
%     if i==1;
%         title([fn sprintf('\n') 'e' num2str(i)]); % allows title to have file name
%         xlabel('frequency (Hz)');
%         ylabel('normalized PSD');
%     elseif i == M1_ch
%         title(['e' num2str(i) ],...
%             'FontWeight','b');
%     else
%         title(['e' num2str(i) ]);
%     end
%     set(ha,'YLim',YLIM_NORM);
%     set(ha,'XLim',XLIM_SPEC1);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
%     
%     % 2nd row subplot
%     
%     subplot(4,5,5+j)
%     ha = gca;
%     hold(ha,'on');
%     plot(freq,normpsdall(:,i),'k','LineWidth',2);
%     if i==1;
%         xlabel('frequency (Hz)');
%         ylabel('normalized PSD');
%     end
%     set(ha,'YLim',YLIM_NORM);
%     set(ha,'XLim',XLIM_SPEC2);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
%     
%     % 3rd row subplot
%     subplot(4,5,2*5+j)
%     ha = gca;
%     hold(ha,'on');
%     plot(freq,logpsdall(:,i),'k','LineWidth',2);
%     if i==1;
%         xlabel('frequency (Hz)');
%         ylabel('log PSD');
%     end
%     set(ha,'YLim',YLIM_LOG);
%     set(ha,'XLim',XLIM_SPEC1);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
%     
%     % 4th row subplot
%     subplot(4,5,3*5+j)
%     ha = gca;
%     hold(ha,'on');
%     plot(freq,logpsdall(:,i),'k','LineWidth',2);
%     if i==1;
%         xlabel('frequency (Hz)');
%         ylabel('log PSD');
%     end
%     set(ha,'YLim',YLIM_LOG);
%     set(ha,'XLim',XLIM_SPEC2);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
% end
% saveas(hf1,[fn '_ecogPSDrest_M1_1'],'fig');
% 
% % third figure plot PSD in 4x6 subplot, 6 columns for ecog/LFP contact around M1, 4 rows for
% % different plot views
% hf1 = figure;
% M1=M1_ch2;
% x=[M1-2 M1-1 M1 M1+1 M1+2];
% for i = x
%     j=i-x(1)+1;
%     % 1st row subplot
%     subplot(4,5,j)
%     ha = gca;
%     hold(ha,'on');
%     plot(freq,normpsdall(:,i),'k','LineWidth',2);
%     if i==1;
%         title([fn sprintf('\n') 'e' num2str(i)]); % allows title to have file name
%         xlabel('frequency (Hz)');
%         ylabel('normalized PSD');
%     elseif i == M1_ch
%         title(['e' num2str(i) ],...
%             'FontWeight','b');
%     else
%         title(['e' num2str(i) ]);
%     end
%     set(ha,'YLim',YLIM_NORM);
%     set(ha,'XLim',XLIM_SPEC1);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
%     
%     % 2nd row subplot
%     
%     subplot(4,5,5+j)
%     ha = gca;
%     hold(ha,'on');
%     plot(freq,normpsdall(:,i),'k','LineWidth',2);
%     if i==1;
%         xlabel('frequency (Hz)');
%         ylabel('normalized PSD');
%     end
%     set(ha,'YLim',YLIM_NORM);
%     set(ha,'XLim',XLIM_SPEC2);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
%     
%     % 3rd row subplot
%     subplot(4,5,2*5+j)
%     ha = gca;
%     hold(ha,'on');
%     plot(freq,logpsdall(:,i),'k','LineWidth',2);
%     if i==1;
%         xlabel('frequency (Hz)');
%         ylabel('log PSD');
%     end
%     set(ha,'YLim',YLIM_LOG);
%     set(ha,'XLim',XLIM_SPEC1);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
%     
%     % 4th row subplot
%     subplot(4,5,3*5+j)
%     ha = gca;
%     hold(ha,'on');
%     plot(freq,logpsdall(:,i),'k','LineWidth',2);
%     if i==1;
%         xlabel('frequency (Hz)');
%         ylabel('log PSD');
%     end
%     set(ha,'YLim',YLIM_LOG);
%     set(ha,'XLim',XLIM_SPEC2);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
% end
% saveas(hf1,[fn '_ecogPSDrest_M1_2'],'fig');
% 
% return;

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
        


