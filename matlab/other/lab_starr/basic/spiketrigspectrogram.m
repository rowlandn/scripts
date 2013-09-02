function spiketrigspectrogram()
% spiketrigecog()
% This m-file takes single-unit spike timestamps stored in a nex, then
% runs spike-triggered signal averaging on ecog recordings.
%
% Created by SAS 7-29-2009

%% Initialize varaibles
Fs = 1000;      % GS4k: sampling rate for ecog
% Fs = 751;      % AlphaOmega: sampling rate for ecog
PRE_SPK = 0.5; % pre-spike trigger period (sec)
PST_SPK = 0.5;  % post-spike trigger period (sec)
STDEV = 3;     % standard deviation for significance
NSPK = 1000;     % number of spikes to use for signal averaging

spk_toggle = false; % true = use NSPK; false = use all spikes in recording

taxis = 1000*(-PRE_SPK:1/Fs:PST_SPK); % time vector, multiplied by 1000 to change to msec scale
faxis = 1:200;  % frequency vector

%% Build a finite impulse response (FIR)bandpass filter kernel with N "taps" 
% and a user-defined bandpass setting.
% (where sampling rate is 1000 Hz, nyquist frequency is
%500 Hz, and the filter bandpass is 0.15 to 0.2 of the distance between 0 and 500).
%the number of elements in the filter kernel is 100 to provide fast
%roll-off.  The filter type is  windowed sinc filter, using the default
%Hamming window

% Define the number of taps
Ntap = 100;    
% Define bandpass filter cutoffs
COlo = 80; % low cutoff of bandpass filter
COhi = 150; % high cutoff of bandpass filter
% Define a 2 element vector with the bandpass window
WN = [COlo,COhi]/(Fs/2);
% Build FIR gamma range bandpass filter
Bg=fir1(Ntap,WN);

% Plot the frequency response for bandpass filter kernel b to show that is
% really is a bandpass filter 
% freqz(B);

%% spike timestamps
[nexfn, nexpn] = uigetfile('*.nex', 'Select spike time file (NEX)');
filename = strrep(nexfn,'.nex','');

if (nexfn ~= 0)
    cd(nexpn);
    [nvar, varname, types] = nex_info(nexfn);
    if nvar > 1      
        % allow user to select desired channel if there are multiple
        % channels
        varnamecell = cellstr(varname);
        varnameidx = menu(['There are > 1 channels in ' nexfn '.' sprintf('\n')...
        'Choose the channel name with desired spike timstamps.'], varnamecell);
        varname = varnamecell{varnameidx};
    end
    [spk.n, spk.t] = nex_ts(nexfn,varname);
else
    error(['I can''t find the NEX file:  ' nexfn ' in ' nexpn]);
end

% find ISI < 1 msec
isi = 1000 .* diff(spk.t);
bad = find(isi < 1);
if ~isempty(bad)
    warning(['Found ' int2str(length(bad)) ' ISI''< 1 msec!!  Correcting...']);
    spk.t( bad+1 ) = [];
    spk.n = length(spk.t);
end

%% ecog data
ecogfn = strcat(filename,'_ecog.mat');
load(ecogfn);

% % Idrissi
% M1_ch=3;

% remontage (units are in microv)
ecog12= (ecog.contact_pair(1).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal);
ecog23= (ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal);
ecog34= (ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal);
ecog45= (ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(5).raw_ecog_signal);
ecog56= ecog.contact_pair(5).raw_ecog_signal;

% pt Solis 3/23/10 recorded with ecog signal and common reference 6 flipped
% around.  correct for flipped signal.
% ecog12= (ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(1).raw_ecog_signal);
% ecog23= (ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal);
% ecog34= (ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal);
% ecog45= (ecog.contact_pair(5).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal);
% ecog56= -ecog.contact_pair(5).raw_ecog_signal;

% % AlphaOmega (SAS 2/9/10): there is an annoying DC offset in all ecog channels.
% % Eliminate offset from remontaged signals by subtracting the average value.
% ecog12=ecog12-mean(ecog12);
% ecog23=ecog23-mean(ecog23);
% ecog34=ecog34-mean(ecog34);
% ecog45=ecog45-mean(ecog45);
% ecog56=ecog56-mean(ecog56);
% 
% % AlphaOmega (SAS 2/9/10): remontaged channels have significantly lower signal amplitude
% % than ecog56, making it difficult to visualize STA waveforms on the plot.
% % Multiply remontaged channels and computationally increase their signal
% % amplitude.
% ecog12=ecog12*10;
% ecog23=ecog23*10;
% ecog34=ecog34*10;
% ecog45=ecog45*10;

ecogall = [ecog12; ecog23; ecog34; ecog45; ecog56];


%% process spike timestamps
spk_t = spk.t;

[ecg_n, ecg_l] = size(ecogall); % find number and length of ecog channels
% discard spk ts at beginning of recording with inadequate pre-spike period
spk_t = spk_t(spk_t > PRE_SPK); 
% discard spk ts at end of recording with inadequate post-spike period
spk_t = spk_t((spk_t+PST_SPK)*Fs < ecg_l); 
% spk_t = spk_t(1:NSPK);
spk_n = length(spk_t);


% optional: use predetermined number of spikes
if spk_toggle
    if spk_n > NSPK
        spk_t(NSPK+1:end) = []; % only include the first NSPK spikes
        spk_n = NSPK;
    else
        error([nexfn ' does not exceed the # spike cutoff of ' num2str(NSPK)...
            '.  Lower minimum cutoff for this unit.']);
    end
end

% %% random timestamps
% tend = spk_t(end)-PST_SPK; % random timestamps are bound in time between 0 sec and time of last spike minus post-spike period
% rnd_t = ceil(1000*(tend*rand(spk_n,1))); % random time (msec)
% rnd_t = rnd_t/1000+PRE_SPK; %random time (sec), plus pre-spike period

%% wavelet transform
% Run KJM's wavelet algorithm on M1 ECoG
tf = dg_tf_pwr(ecogall(M1_ch,:));

%% M1 bandpass, hilbert transform
M1bpg=conv(ecogall(M1_ch,:),Bg,'same');
M1bpgaa=abs(hilbert(M1bpg));
% M1bpgaax2=abs(hilbert(M1bpgaa));

%% time-reversed signals
% time-reverse tf and BP signal for calculating confidence interval, using only
% the recording segment from which the spikes were taken
i_end = int32((spk_t(end)+PST_SPK)*Fs);
tff=flipud(tf(1:i_end,:));
M1bpgaaf=fliplr(M1bpgaa(1:i_end));

%% process spike-triggerd data

% initialize
A = zeros(length(taxis),length(faxis)); %2D matrix for storing wavelet transformed data
Af = zeros(length(taxis),length(faxis)); %2D matrix for storing time-reversed wavelet transformed data
B = zeros(1,length(taxis)); %1D matrix for storing gamma analytic amplitude data
Bf = zeros(1,length(taxis)); %1D matrix for storing time-reversed gamma analytic amplitude data

% B2 = zeros(1,length(taxis)); %1D matrix for storing gamma analytic amplitude x2 data

% parse each stim epoch
for i = 1:spk_n
    tmp1 = int32((spk_t(i)-PRE_SPK) * Fs);  % int32 used to keep index in integer format
    tmp2 = int32((spk_t(i)+PST_SPK) * Fs);
    % since int32 rounds to the next closest integer, there may be some
    % cases where the parsing indeces are off by 1.  fix them on the fly.
    d = tmp2 - tmp1;
    % lengt(t)-d must equal 1 for the parsed ecog data to fit S 3D matrix
    if length(taxis)-d == 0
        tmp2 = tmp2-1;
    elseif length(taxis)-d == 2
        tmp2 = tmp2+1;
    end
    A = A + tf(tmp1:tmp2,:);
    Af = Af + tff(tmp1:tmp2,:);
    B = B + M1bpgaa(tmp1:tmp2);
    Bf = Bf + M1bpgaaf(tmp1:tmp2);
%     B2 = B2 + M1bpgaax2(tmp1:tmp2);
end

A=A/spk_n;
Af=Af/spk_n;
B=B/spk_n;
Bf=Bf/spk_n;
% B2=B2/spk_n;

%% process random-triggered data
% 
% Ar = zeros(length(taxis),length(faxis)); %2D matrix for storing wavelet transformed data
% Br = zeros(1,length(taxis)); %1D matrix for storing gamma analytic amplitude data
% % Br2 = zeros(1,length(taxis)); %1D matrix for storing gamma analytic amplitude x2 data
% % parse each randomly triggered epoch
% for i = 1:spk_n
%     tmp1 = int32((rnd_t(i)-PRE_SPK) * Fs);  % int32 used to keep index in integer format
%     tmp2 = int32((rnd_t(i)+PST_SPK) * Fs);
%     % since int32 rounds to the next closest integer, there may be some
%     % cases where the parsing indeces are off by 1.  fix them on the fly.
%     d = tmp2 - tmp1;
%     % lengt(t)-d must equal 1 for the parsed ecog data to fit S 3D matrix
%     if length(taxis)-d == 0
%         tmp2 = tmp2-1;
%     elseif length(taxis)-d == 2
%         tmp2 = tmp2+1;
%     end
%     Ar = Ar + tf(tmp1:tmp2,:);
%     Br = Br + M1bpgaa(tmp1:tmp2);
%     Br2 = Br2 + M1bpgaax2(tmp1:tmp2);
% end
% 
% Ar=Ar/spk_n;
% Br=Br/spk_n;
% % Br2=Br2/spk_n;

%% modulation index
mi1 = (mean(B(301:700))-mean(Bf(301:700)))/std(Bf(301:700));
% also calcuate modulation index using DC corrected, rectified signal
B2 = abs(B-mean(B));
Bf2 = abs(Bf-mean(Bf));
% mi2 = (mean(B2(401:600))-mean(Bf2(401:600)))/std(Bf2(401:600));
mi2 = (mean(B2(401:500))-mean(Bf2(401:500)))/std(Bf2(401:500));

%% plot figures
% spike-triggered M1 spectrogram
hf1=figure;
valmin=min(min(A));
valmax=max(max(A));
clims=[valmin valmax];
subplot(2,1,1);
imagesc(taxis,faxis,A',clims);
hold on;
plot([0 0],[1 200],'k--');
title([filename sprintf('\n')...
    'M1 LFP spike triggered, #spks=' num2str(spk_n)],...
    'FontSize',12);
% xlabel('Time (msec)');
ylabel('Frequency (Hz)');

% spike-triggered time-reversed M1 spectrogram
% hf2=figure;
subplot(2,1,2);
imagesc(taxis,faxis,Af',clims);
hold on;
% title([filename sprintf('\n')...
%     'M1 LFP time-reversed, #spks=' num2str(spk_n)],...
%     'FontSize',12);
plot([0 0],[1 200],'k--');
title('M1 LFP time-reversed',...
    'FontSize',12);
xlabel('Time (msec)');
ylabel('Frequency (Hz)');

colorbar([0.9307 0.1048 0.02354 0.8226]);

% hf3=figure;
% % create a vector of ones
% vones = ones(1,length(taxis));
% Bmean=vones*mean(B(401:600));
% Bfmean=vones*mean(Bf(401:600));
% thu = Bfmean + (STDEV * std(Bf(401:600)));
% thl = Bfmean - (STDEV * std(Bf(401:600)));
% plot(taxis,B,'k','LineWidth',2),hold on,plot(taxis,Bf);
% plot(taxis,Bmean,'-.k','LineWidth',2);
% plot(taxis,Bfmean,'-.');
% plot(taxis,thu,'b-',taxis,thl,'b-');
% ylm=ylim;
% plot([0 0],[ylm(1) ylm(2)],'k--');
% xlim([-200 200]);
% % text(100,thu(1),...
% %     ['modulation index 1= ' num2str(mi1)],...
% %     'HorizontalAlignment','left',...
% %     'VerticalAlignment','bottom');
% title([filename sprintf('\n')...
%     'M1 gamma BP AA, #spks=' num2str(spk_n) ',sig thresh=' num2str(STDEV) 'SD'  sprintf('\n')...
%     'modulation index=' num2str(mi1) ', -100 to 0msec'],...
%     'FontSize',12);
% legend('spike-triggered','time-reversed');
% legend('boxoff');

hf4=figure;
% create a vector of ones
vones = ones(1,length(taxis));
B2mean=vones*mean(B2(401:500));
Bf2mean=vones*mean(Bf2(401:500));
% create confidence intervals
% th2=Bf2mean+(2*std(Bf2(401:500)));
th3=Bf2mean+(3*std(Bf2(401:500)));
th4=Bf2mean+(4*std(Bf2(401:500)));
th5=Bf2mean+(5*std(Bf2(401:500)));

plot(taxis,B2,'k','LineWidth',2),hold on,plot(taxis,Bf2);
plot(taxis,B2mean,'-.k','LineWidth',2);
plot(taxis,Bf2mean,'-.');
% plot(taxis,th2,'b-');
plot(taxis,th3,'b-');
plot(taxis,th4,'b-');
plot(taxis,th5,'b-');
ylm=ylim;
plot([0 0],[ylm(1) ylm(2)],'k--');
xlim([-200 200]);
% text(100,th(1),...
%     ['modulation index 2= ' num2str(mi2)],...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','bottom');
% text(180,th2(1),...
%     ['2SD'],...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','bottom');
text(180,th3(1),...
    ['3SD'],...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(180,th4(1),...
    ['4SD'],...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(180,th5(1),...
    ['5SD'],...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
title([filename sprintf('\n')...
    'M1 gamma BP AA, #spks=' num2str(spk_n) sprintf('\n')...
    'modulation index=' num2str(mi2) ' (-100 to 0msec)' ],...
    'FontSize',12);
legend('spike-triggered','time-reversed');
legend('boxoff');

% 1-30Hz frequency view
% figure,imagesc(taxis,faxis(1:30),E(:,1:30)');
% colorbar;
% title([filename sprintf('\n')...
%     'M1 LFP, #spks=' num2str(spk_n)],...
%     'FontSize',12);
% xlabel('Time (msec)');
% ylabel('Frequency (Hz)');

%% save variables and figures
% saveas(hf1,[filename '_wvlt'],'fig');
% save([filename '_wvlt'],'taxis','faxis','tf','E','spk_t','spk_n');

