function ecogfft5
% this function plots the frequency spectrum of files containing digitized
% ecog data.

% define variables
FREQ_LO = [8 30];   % beta band
FREQ_MED = [35 57];  %low gamma band
FREQ_HI = [78 100]; % high gamma band
Fs=1000;            % samples per unit time, in this case digitization at 1000 hz rate
OFFSET = 0.5;        % epoch time offset in seconds
EPOCH_LEN = 2;       % epoch length in seconds 
SPEC_MAX = 100;     % max end of spectral frequency for plotting
SPAN = 9;          % span of convolution kernel
YLIM = [-50 50];    % voltage level to be graphed for raw ecog signal 

%%
% import filename and pathname of mat file containing ecog data created by APMconv5.m
[fn pn] = uigetfile('*.mat','Select .mat containing ecog data');
cd(pn);
% load ecog_rest and ecog_active variables and store both into one strucutre called ecog
load([pn fn]);
% remove '_ecog.mat' ending from filename
fn = strrep(fn,'_ecog.mat','');

hf = figure;
num_contact_pair = length(ecog.contact_pair);
num_epoch = length(ecog.rest_time);

% initialize structures that will contain all rest/active analysis
rest = struct('contact_pair',{});
active = struct('contact_pair',{});

           
for i = 1:num_contact_pair
    
    % parse each rest/active epoch from each contact pair
    
    for j = 1:num_epoch
        start_rest = (ecog.rest_time(j) + OFFSET) * Fs; % time offset added to epoch times
        end_rest = start_rest + (Fs * EPOCH_LEN) - 1;
        rest(1).contact_pair(i).epoch(j).raw_ecog_signal = ...
            ecog.contact_pair(i).raw_ecog_signal(start_rest:end_rest);
        start_active = (ecog.active_time(j)+ OFFSET) * Fs;
        end_active = start_active + Fs*EPOCH_LEN - 1;
        active(1).contact_pair(i).epoch(j).raw_ecog_signal = ...
            ecog.contact_pair(i).raw_ecog_signal(start_active:end_active);
    end
    
%     %graph raw signal, both rest/active from epoch 1 ONLY
%     
%     t_raw_ecog = 1: Fs*EPOCH_LEN;
%     subplot(6,num_contact_pair,i); % 1st row subplots
%     plot(t_raw_ecog, rest.contact_pair(i).epoch(1).raw_ecog_signal);
%     ylim(YLIM);
%     if i == 1
%         ylabel(['raw ecog signal at rest' sprintf('\n')...
%             'epoch 1'],...
%             'FontSize',8);
%         xlabel('time (ms)');
%     end
%     if i <= num_contact_pair
%         GL4K_channel_number = ['Channel #'...
%             num2str(ecog.contact_pair(i).chan_num)];
%         title(GL4K_channel_number);
%     end
%     
%     subplot(6,num_contact_pair,num_contact_pair+i);  % 2nd row subplots
%     plot(t_raw_ecog, active.contact_pair(i).epoch(1).raw_ecog_signal);
%     ylim(YLIM);
%     if i == 1
%         ylabel(['raw ecog signal' sprintf('\n')...
%             'during active movement' sprintf('\n')...
%             'epoch 1'],...
%             'FontSize',8);
%         xlabel('time (ms)');
%     end
    
% graph averaged unsmoothed FFT, superimosed rest/active
    for j = 1:num_epoch  %perform fft for each epoch
        [rest(1).contact_pair(i).epoch(j).unsmoothed_norm_power_fft ...
            rest(1).contact_pair(i).epoch(j).raw_power_fft] =...
            fft_ecog_signal(rest.contact_pair(i).epoch(j).raw_ecog_signal);
        [active(1).contact_pair(i).epoch(j).unsmoothed_norm_power_fft ...
            active(1).contact_pair(i).epoch(j).raw_power_fft] =...
            fft_ecog_signal(active.contact_pair(i).epoch(j).raw_ecog_signal);
    end

   % calculate mean across all epochs
   % note: mean of 'unsmoothed_norm_power_fft' is used for graphing whereas
   % mean of 'raw_power_fft' is for data anlysis 
    rest(1).contact_pair(i).mean_unsmoothed_norm_power_fft = ...
        mean_signal(rest.contact_pair(i),'unsmoothed_norm_power_fft'); % calculate mean for fieldname 'unsmoothed_norm_power_fft'
    active(1).contact_pair(i).mean_unsmoothed_norm_power_fft = ...
        mean_signal(active.contact_pair(i),'unsmoothed_norm_power_fft');
    rest(1).contact_pair(i).mean_raw_power_fft = ...
        mean_signal(rest.contact_pair(i),'raw_power_fft'); % calculate mean for fieldname 'raw_power_fft'
    active(1).contact_pair(i).mean_raw_power_fft = ...
        mean_signal(active.contact_pair(i),'raw_power_fft');
    
    freq_unsmoothed = linspace(0,Fs/2,Fs*EPOCH_LEN/2);
    
    subplot(4,num_contact_pair,i); % 1st row subplots
    ha = gca;       % create handle for current axes
    hold(ha,'on');
    plot(freq_unsmoothed, rest.contact_pair(i).mean_unsmoothed_norm_power_fft,...
        freq_unsmoothed, active.contact_pair(i).mean_unsmoothed_norm_power_fft,'r');
    GL4k_chan_num = num2str(ecog.contact_pair(i).chan_num);
    if i <= num_contact_pair && i~=1
        title(['Channel #' GL4k_chan_num]);
    elseif i == 1
        title([fn sprintf('\n') 'Channel #' GL4k_chan_num]);
        ylabel(['average of normalized power' sprintf('\n')...
            'across all epochs'],...
        'FontSize',8);
        xlabel('frequency (Hz)');
    end
    hl = legend('Rest','Active');
    set(hl,'FontSize',5,'Box','off');
    set(ha,'YLim',[0 1]);
    xlim([0 SPEC_MAX]);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
   
    % graph averaged smoothed FFT, superimposed rest/active
    
    for j = 1:num_epoch     % smooth each epoch using conv_norm_power
        rest.contact_pair(i).epoch(j).convn_norm_power_fft = ...
            convn_power(rest.contact_pair(i).epoch(j).unsmoothed_norm_power_fft, freq_unsmoothed,SPEC_MAX,SPAN);
        active.contact_pair(i).epoch(j).convn_norm_power_fft =...
            convn_power(active.contact_pair(i).epoch(j).unsmoothed_norm_power_fft, freq_unsmoothed, SPEC_MAX,SPAN);
    end
    
    rest(1).contact_pair(i).mean_convn_norm_power_fft = ...
        mean_signal(rest.contact_pair(i),'convn_norm_power_fft');
    active(1).contact_pair(i).mean_convn_norm_power_fft = ...
        mean_signal(active.contact_pair(i),'convn_norm_power_fft');
    
    freq_convn = 0:SPEC_MAX-1;
    
    subplot(4,num_contact_pair,num_contact_pair+i); % 2nd row subplots
    ha = gca;
    hold(ha,'on');
    plot(freq_convn,rest.contact_pair(i).mean_convn_norm_power_fft,...
        freq_convn,active.contact_pair(i).mean_convn_norm_power_fft,...
        'r','LineWidth',2);
    if i == 1
        ylabel(['average of smoothed' sprintf('\n')...
            'normalized power' sprintf('\n')...
            'across all epochs'],...
            'FontSize',8);
        xlabel('frequency (Hz)');
    end
    hl = legend('Rest','Active');
    set(hl,'FontSize',5,'Box','off');
    set(ha,'YLim',[0 1]);
    xlim([0 SPEC_MAX]);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');

    
    % graph log averaged smoothed FFT, superimposed rest/active
    subplot(4,num_contact_pair,2*num_contact_pair+i); % 3rd row subplot
    ha = gca;
    hold(ha,'on');
    rest.contact_pair(i).log_mean_convn_norm_power_fft = ...
        log(rest.contact_pair(i).mean_convn_norm_power_fft);
    active.contact_pair(i).log_mean_convn_norm_power_fft = ...
        log(active.contact_pair(i).mean_convn_norm_power_fft);
    plot(freq_convn, rest.contact_pair(i).log_mean_convn_norm_power_fft,...
        freq_convn, active.contact_pair(i).log_mean_convn_norm_power_fft,...
        'r','LineWidth',2);
    if i == 1
        ylabel(['log average of smoothed' sprintf('\n')...
            'normalized power' sprintf('\n')...
            'across all epochs'],...
            'FontSize',8);
        xlabel('frequency (Hz)');
    end
    hl = legend('Rest','Active');
    set(hl,'FontSize',5,'Box','off');
    xlim([0 SPEC_MAX]);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma bands
    hold(ha,'off');
    
    %graph normalized power

    ratio = 100* (active.contact_pair(i).mean_convn_norm_power_fft./rest.contact_pair(i).mean_convn_norm_power_fft);
    subplot(4,num_contact_pair,3*num_contact_pair+i); % 4th row subplot
    ha = gca;
    hold(ha,'on');
    plot(freq_convn, ratio,'r','LineWidth',2)
    plot([0 SPEC_MAX],[100 100],'--k');
    if i == 1
        ylabel('normalized power (%)',...
            'FontSize',8);
        xlabel('frequency (Hz)');
    end
    hl = legend('Active');
    set(hl,'FontSize',5,'Box','off');
    xlim([0 SPEC_MAX]);
    if max(ratio)<1000  % cut off y-axis limit at 1000
        ylim([0 max(ratio)+50]);
    else
        ylim([0 1000]);
    end
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma bands
    hold(ha,'off');
    
%     % graph normalized power (using non-normalized raw power output of fft)
%     a = convn_power(active.contact_pair(i).mean_raw_power_fft,freq_unsmoothed,SPEC_MAX,SPAN);
%     b = convn_power(rest.contact_pair(i).mean_raw_power_fft,freq_unsmoothed,SPEC_MAX,SPAN);
%     ratio = 100*(a./b);
%     
%     %     ratio = 100* (active.contact_pair(i).mean_raw_power_fft./rest.contact_pair(i).mean_raw_power_fft);
%     %     ratio = convn_power(ratio,freq_unsmoothed,SPEC_MAX,SPAN);
%     
%     subplot(5,num_contact_pair,4*num_contact_pair+i); % 5th row subplot
%     ha = gca;
%     hold(ha,'on');
%     plot(freq_convn, ratio,'r','LineWidth',2)
%     plot([0 SPEC_MAX],[100 100],'--k');
%     if i == 1
%         ylabel('normalized power (%)',...
%             'FontSize',8);
%         xlabel('frequency (Hz)');
%     end
%     hl = legend('Active Movement');
%     set(hl,'FontSize',5,'Box','off');
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma bands
%     xlim([0 SPEC_MAX]);
%     ylim([0 200]);
%     hold(ha,'off');
    
end



% % calculate total power output for rest/active in each frequency band
% FREQ_ALL = [FREQ_LO; FREQ_MED; FREQ_HI];
% active_rest_ratio = [];
% for i = 1:num_contact_pair
%     for j = 1:length(FREQ_ALL)
%         tmp = rest.contact_pair(i).log_mean_convn_norm_power_fft(FREQ_ALL(j,1):FREQ_ALL(j,2));
%         rest(1).contact_pair(i).total_power(j) = sum(tmp);
%         tmp = active.contact_pair(i).log_mean_convn_norm_power_fft(FREQ_ALL(j,1):FREQ_ALL(j,2));
%         active(1).contact_pair(i).total_power(j) = sum(tmp);
%         active_rest_ratio(j,i) = active(1).contact_pair(i).total_power(j)/rest(1).contact_pair(i).total_power(j);
%     end
% end

% save figure
saveas(hf,[fn '_ecog'],'fig');


function [fft_out P_out] = fft_ecog_signal(raw_ecog_signal)
% fft_ecog_signal performs FFT on raw ecog signal and outputs normalized power
% in the frequency domain up to the Nyquist frequency

Y=fft(raw_ecog_signal);
% will use syntax Y=fft(ecogsnippet1, N) where N is the number of points that are
% used for the time domain (and frequency domain), use N=2000)
% note in this section we will take the starting times for each rest
% snippet or each active snippet, run the fft for each.
Y(1)=0; %for some reason the first number in the fft vector is the sum of all others,
% and to avoid confusion in the plotting, we set that first number to zero
n=length(Y);
P=(abs(Y).^2)/n; %calculates the power of the fft, which is square of its magnitude
P_out = P(1:floor(n/2));
% at this point, do the quantitative analyses - many of the functions below
% this point (normalizing to peak power, smoothing, taking logs) are really
% just to facilitate plotting and data visualization

%lets try the following for each individual channel, in each individual movement type (jaw, hand etc)
%The five "raw" P vectors for each epoch should be averaged together (prior to
%any normalization).  then "normalize" the movement related averaged P vector by the averaged rest
%P vector. Plot normalized P vector analagous to fig 2 of liu and aziz 2008) Then segment the normalized P vector into frequency bins
%(bandwidth to be determined, but could start with 10 hz bands 0 to 150)
%and determine mean normalized power in that frequency bin.  write these
%data to a table as in the example.  
%rest vector 

norm_P=P./max(P); %normalizes the power to the maximum power of 1 to facilate plotting
fft_out = norm_P(1:floor(n/2));%pulls out only the first half of the spectrum up to nyquist frequency
% n1 = length(fft_out);
% freq_out =(0:n1-1)*(Fs/n); %sets the frequency or x scale to be correct range to correspond to fft
% % run the above sequence for each of the 5-10 snippets at rest, and each in
% % movement, and add the fft_out vectors for each to a final
return

function out = mean_signal(contact_pair, fieldname)
num_epoch = length(contact_pair.epoch);
A = [];
for j = 1:num_epoch
    tmp = eval(['contact_pair.epoch(' num2str(j) ').' fieldname]);
    A = [A;tmp]; % concatenate
end
out = mean(A);
return


function out = convn_power(power_fft,frequency,SPEC_MAX,SPAN)

edge = zeros(1,floor(SPAN/2));
out = [edge power_fft edge];  % add edges to avoid edge effect during convolution
k = (1/SPAN)*ones(1,SPAN);  % create convolution kernel

out = convn(out,k); % convolve
% remove edges and extra elements added by convn
% note: convn adds k-1 elements to output array, see doc convn for details
out = out(2*length(edge)+1:length(out)-2*length(edge));
% downsample, only sample at non-rational (whole number) frequency
frequency = floor(frequency); % sample at non-rational number frequency
frequency = frequency(frequency<SPEC_MAX); % extract the frequencies below SPEC_MAX

idx = 1;
for i=1:length(frequency)-1
    if frequency(i)+1 == frequency(i+1)
        idx = [idx i+1];
    end
end
out = out(idx);
return

% function out = downsample_fft(power_fft,frequency,SPEC_MAX)
% % downsample, only sample at non-rational (whole number) frequency
% frequency = floor(frequency); % sample at non-rational number frequency
% frequency = frequency(frequency<SPEC_MAX); % extract the frequencies below SPEC_MAX
% 
% idx = 1;
% for i=1:length(frequency)-1
%     if frequency(i)+1 == frequency(i+1)
%         idx = [idx i+1];
%     end
% end
% out = power_fft(idx);
% return

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
return




