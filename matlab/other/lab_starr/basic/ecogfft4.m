function ecogfft4
% this function plots the frequency spectrum of files containing digitized
% ecog data.


% define variables
FREQ_LO = [8 30];   % beta band
FREQ_MED = [35 57];  %low gamma band
FREQ_HI = [78 100]; % high gamma band
Fs=1000;            % samples per unit time, in this case digitization at 1000 hz rate
BIN = 1;            % bin used for FFT downsampling
SPEC_MAX = 100;     % max end of spectral frequency for plotting
SPAN = 9;          % span of convolution kernel
YLIM = [-50 50];  % for graphing

%%
% import filename and pathname of mat file containing ecog rest data created by APMconv5.m
[fn_rest pn_rest] = uigetfile('*.mat','First select mat file containing ecog at rest');
cd(pn_rest);
    % standardize signal length in time (2 sec), possibly take snippets of
    % recordings in the future
% import filename and pathname of mat file containing ecog active data
[fn_active pn_active] = uigetfile('*.mat','Now select mat file containing ecog during active movement');

% load ecog_rest and ecog_active variables and store both into one strucutre called ecog
load([pn_rest fn_rest]);
ecog.rest.contact_pair = ecog_rest(1).contact_pair;
load([pn_active fn_active]);
ecog.active.contact_pair = ecog_active(1).contact_pair;
% may 12 - will change to recording a single apm file for each type of
% movement (eg rest/active hand will be one apm, rest/active elbow for
% another apm, rest/active shoulder for another 
% each apm will have multiple epochs of rest and movement.  the emg or
% accel channel will show when those epochs start and end.  starr will
% examine the files first and determine the time (or sample #) of the start
% of each epoch, eg rest periods start at samples 0, 5000, 9000, ,
% while active start at samples 3000, 7500, 1200, etc.
% for graphing the time domain signal, you can just graph the first time
% snippets for rest and active (for each channel) rather than graph the
% entire time domain file.

% create for-loop that will graph each contact pair
hf = figure;
num_contact_pair = length(ecog.rest.contact_pair);
for i = 1:num_contact_pair
    
    %graph raw signal, both rest/active
    raw_ecog_signal_rest = ecog.rest.contact_pair(i).raw_ecog_signal;
    z =length(raw_ecog_signal_rest);
    lastnum=(z-1)/Fs;
    t_ecog_rest=0:1/Fs:lastnum;
    subplot(5,num_contact_pair,i);
    plot(t_ecog_rest, raw_ecog_signal_rest);
    ylim(YLIM);
    if i == 1
        ylabel('raw ecog signal at rest',...
            'FontSize',8);
    end
    if i <= num_contact_pair
        GL4K_channel_number = ['Channel #' num2str(i+1)];
        title(GL4K_channel_number);
    end
    
    raw_ecog_signal_active = ecog.active.contact_pair(i).raw_ecog_signal;
    z =length(raw_ecog_signal_active);
    lastnum=(z-1)/Fs;
    t_ecog_active=0:1/Fs:lastnum;
    subplot(5,num_contact_pair,num_contact_pair+i);
    plot(t_ecog_active, raw_ecog_signal_active);
    ylim(YLIM);
    if i == 1
        ylabel(['raw ecog signal',sprintf('\n'), 'during active movement'],...
            'FontSize',8);
    end

    
    % graph unsmoothed FFT, superimosed rest/active
    % sumperimposing is currently not possible because unsmoothed rest/active fft lengths do not
    % match, this will be fixed once signal length is standardized
    [ecog.rest.contact_pair(i).unsmoothed_fft freq_rest]=...
        fft_ecog_signal(ecog.rest.contact_pair(i).raw_ecog_signal,Fs);
    [ecog.active.contact_pair(i).unsmoothed_fft freq_active]=...
        fft_ecog_signal(ecog.active.contact_pair(i).raw_ecog_signal,Fs);
    
    % graph smoothed FFT, superimposed rest/active
    ecog.rest.contact_pair(i).smoothed_fft = ...
        smooth_norm_power(ecog.rest.contact_pair(i).unsmoothed_fft, freq_rest,SPEC_MAX,BIN);
    ecog.active.contact_pair(i).smoothed_fft =...
        smooth_norm_power(ecog.active.contact_pair(i).unsmoothed_fft, freq_active, SPEC_MAX,BIN);
    
    % optionally, use convn_norm_power to convolve then downsample
%     ecog.rest.contact_pair(i).smoothed_fft = ...
%         convn_norm_power(ecog.rest.contact_pair(i).unsmoothed_fft, freq_rest,SPEC_MAX,BIN,SPAN);
%     ecog.active.contact_pair(i).smoothed_fft =...
%         convn_norm_power(ecog.active.contact_pair(i).unsmoothed_fft, freq_active, SPEC_MAX,BIN,SPAN);
    

    freq_smooth = 0:BIN:SPEC_MAX-BIN;  % frequency domain after fft normalized power has been smoothed
    subplot(5,num_contact_pair,3*num_contact_pair+i);
    ha1 = gca;
    hold(ha1,'on');
    plot(freq_smooth,ecog.rest.contact_pair(i).smoothed_fft,...
        freq_smooth,ecog.active.contact_pair(i).smoothed_fft,'LineWidth',2);
    legend1 = legend('Rest','Active Movement');
    set(legend1,'FontSize',5,'Box','off');
    if i == 1
        ylabel('Smoothed FFT',...
            'FontSize',8);
    end
    FFT_band_fill(ha1,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma bands
    hold(ha1,'off');
    
    % graph log smoothed FFT, superimposed rest/active
    subplot(5,num_contact_pair,4*num_contact_pair+i);
    ha2=gca;
    hold(ha2,'on');
    plot(freq_smooth,log(ecog.rest.contact_pair(i).smoothed_fft),...
        freq_smooth,log(ecog.active.contact_pair(i).smoothed_fft),'LineWidth',2);
    legend2 = legend('Rest','Active Movement');
    set(legend2,'FontSize',5,'Box','off');
    if i == 1
        ylabel('Log smoothed FFT',...
            'FontSize',8);
    end
    FFT_band_fill(ha2,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma bands
    hold(ha2,'off');
    
end


function [fft_out freq_out] = fft_ecog_signal(raw_ecog_signal,Fs)
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
P=P./max(P); %normalizes the power to the maximum power of 1 to facilate plotting
fft_out = P(1:floor(n/2));%pulls out only the first half of the spectrum up to nyquist frequency
n1 = length(fft_out);
freq_out =(0:n1-1)*(Fs/n); %sets the frequency or x scale to be correct range to correspond to fft
% run the above sequence for each of the 5-10 snippets at rest, and each in
% movement, and add the fft_out vectors for each to a final
return

function out = smooth_norm_power(norm_power, frequency,SPEC_MAX,BIN)
% smooth_norm_power smoothes normalized power data

% initialize
out = zeros(1,SPEC_MAX/BIN);

for i = 1:length(out)
    I_bin = (frequency>=BIN*(i-1)) & (frequency<BIN*i); % the bins will be left-edge inclusive
    out(i) = mean(norm_power(I_bin));   % smooth by calculating average values inside each bin,
    % maybe change to convolution and downsample in the future
end
return

function out = convn_norm_power(norm_power,frequency,SPEC_MAX,BIN,SPAN)

edge = zeros(1,floor(SPAN/2));
out = [edge norm_power edge];  % add edges to avoid edge effect during convolution
k = (1/SPAN)*ones(1,SPAN);  % create convolution kernel

out = convn(out,k); % convolve
% remove edges and extra elements added by cnvn
% note:cnvn adds k-1 elements to output array, see doc cnvn for details
out = out(2*length(edge)+1:length(out)-2*length(edge));
% downsample
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

function FFT_band_fill(ha,FREQ_HI,FREQ_LO)
% FFT_band_fill fills bands of specified frequncies on axes specified by
% the handle ha
set(ha,'XLim',[0 100]);
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




