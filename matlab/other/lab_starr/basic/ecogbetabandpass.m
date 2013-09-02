function ecogbetabandpass
%this program will bandpass filter an ecog signal with a 4 hz band centered
%on the peak frequency

Fs=1000;                % samples per unit time, in this case digitization at 1000 hz rate
% Fs = 1.502403855323792e3;   % Alpha Omega ecog/lfp are recorded using different sampling rate. 
% note: in the future, consider downsampling AO data to make it as close to
% GL4k data as possible
OFFSET = 0.5;           % epoch time offset in seconds
EPOCH_LEN = 2;          % epoch length in seconds (1.792 is the lenthg of 6 overlapping segments)
WINDOW = 512;           % segment length and Hamming window length for welch's method
NOVERLAP = 256;         % # signal samples that are common to adjacent segments for welch's method
NFFT = 512;             % length of fft
[fn pn] = uigetfile('*_ecog.mat','Select _ecog.mat containing ecog/lfp data');
cd(pn);
% load ecog data
load([pn fn]);
%define variable m1ecog which is the desired signal to analyze

% Determine GS4000 or AlphaOmega, set correction factor

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
ecog12=ecog.contact_pair(1).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal;
ecog23=ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal;
ecog34=ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal;
ecog45=ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(5).raw_ecog_signal;
ecog56=ecog.contact_pair(5).raw_ecog_signal;
 
ecogall = [ecog12; ecog23; ecog34; ecog45; ecog56];
if length(ecog.contact_pair)==6
    % assume LFP data is in contact pair 6
    lfp=ecog.contact_pair(6).raw_ecog_signal;
    ecogall = [ecogall; lfp];
end

necog = size(ecogall,1);

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
    % correct for frequency response
    for j=1:length(crxnfactor)
        psdall(j,i)=psdall(j,i)/(crxnfactor(j)^2); 
    end
end
  
if ~isempty(M1_ch)
    logpsdall = log10(psdall);
    norm_idx = find(freq>8 & freq<100);
    beta_idx = find(freq>13 & freq<30);
    idx = [];
    
    for i = 1:necog
        if i == M1_ch
            idx = [idx i];
        end
    end
          
        
        maxpower = max(psdall(beta_idx(1):beta_idx(end), idx(1))); 
        %maxfreq = [freq(power_idx1) freq(power_idx2)];
        maxfreq = freq(psdall(:,idx(1)) == maxpower(1));
        
end

peakfreq = maxfreq(1); %modify this code to define peak frequency based on max of the PSD of the ecog signal
lowcut = (peakfreq-2)./500;
highcut = (peakfreq+2)./500;

%Define a finite impulse response (FIR)bandpass filter kernel with 100 "taps" 
% and a bandpass of (peakfreq-2)Hz to (peakfreq + 2 Hz)
% (where sampling rate is 1000 Hz, nyquist frequency is
%500 Hz, and the filter bandpass is lowcut to highcut of the distance between 0 and 500).
%the number of elements in the filter kernel is 100 to provide fast
%roll-off.  The filter type is  windowed sinc filter, using the default
%Hamming window

%define a 2 element vector with the bandpass window
wn=[lowcut,highcut];
b=fir1(100, wn); %defines a bandpass filter kernal b that will perform the desired bandpass filtering when implemented by convolution

%NOTE: WE MAY NEED TO FILTER BIRECTIONNALLY TO REMOVE PHASE SHIFTS, OR USE
%EEGFILT
%plot the frequency response for bandpass filter kernel b to show that is
%really is a bandpass filter 
freqz(b);



%implement the FIR filter kernel b described above, by convolution. 
% the convolution "shape" is "same", which returns a signal the same length
% as the original time signal s
%plot the original signal and its beta bandpassed filtered version

%M1_ch_beta=conv(ecogall(M1_ch,:), b); %Might take a long time on a long file
M1_ch_beta = filtfilt(b,1,ecogall(M1_ch,:));
plot(M1_ch_beta);
axis([0 1000 -400 400]); %plots the first 1 sec of data from -400 to +400 microvolts
hold on;
plot(ecogall(M1_ch,:), 'r');


%Find the beta phase from the angle of the complex analytic function for M1ecogbeta

aampbeta=(abs(hilbert(M1_ch_beta)));
aphasebeta=(angle(hilbert(M1_ch_beta)));
tr_time_idx = find(aphasebeta < -3);
idx1 = diff(idx)<30;
tr_time_idx(idx1) = [];  

T = 1/Fs;
t = 0:T:T*(length(aphasebeta)-1);
tr_time = t(tr_time_idx);
        
        fn = strrep(fn,'_ecog.mat','');
        save([fn, '_beta_trough_data.mat'],'ecog','ecogall','tr_time', 'M1_ch','peakfreq', 'M1_ch_beta', 'aampbeta', 'aphasebeta');
end

% find the times of the troughs of the beta oscillations 
%NEXT USE THE TIMESTAMPS OF BETA TROUGHS IN THE WAVELET SPECTROGRAM
%PROGRAM.


