function spiketrigecog2()
% spiketrigecog()
% This m-file takes single-unit spike timestamps stored in a nex, then
% runs spike-triggered signal averaging on ecog recordings.
%
% Created by SAS 7-29-2009

% spiketrigecog2 uses a small window around t=0 to calculate standard
% deviation.

%% Initialize varaibles
Fs = 1000;      % GS4k: sampling rate for ecog
% Fs = 751;      % AlphaOmega: sampling rate for ecog
PRE_SPK = 0.5; % pre-spike trigger period (sec)
PST_SPK = 0.5;  % post-spike trigger period (sec)
STDEV = 2.5;     % standard deviation for significance
NSPK = 2000;     % number of spikes to use for signal averaging
WINDOW = 512;           % segment length and Hamming window length for welch's method
NOVERLAP = 256;         % # signal samples that are common to adjacent segments for welch's method
NFFT = 512;             % length of fft

spk_toggle = true; % true = use NSPK; false = use all spikes in recording


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

% load data
% [ecogfn, ecogpn] = uigetfile('*_ecog.mat', 'Select ecog file (*_ECOG.MAT)');
% cd(ecogpn);
ecogfn = strcat(filename,'_ecog.mat');
load(ecogfn);

% remontage (units are in microv)
ecog12= (ecog.contact_pair(1).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal);
ecog23= (ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal);
ecog34= (ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal);
ecog45= (ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(5).raw_ecog_signal);
ecog56= ecog.contact_pair(5).raw_ecog_signal;

M1_ch=3;
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


%% spike-triggered average of ECog signal

spk_t = spk.t;

% if clause below for using burst onset to analyze units with tremor
% frequency band activity
% if false
%     Burst = legendy_new3(isi,[],[],[],[],5);
%     spk_t =ones(1,length(Burst));
%     for i=1:length(Burst)
%         spk_t(i) = sum(isi(1:Burst(i).begin));
%     end
%     spk_t = spk_t/1000; % convert from msec -> sec
% end

[ecg_n, ecg_l] = size(ecogall); % find number and length of ecog channels
% discard spk ts at beginning of recording with inadequate pre-spike period
spk_t = spk_t(spk_t > PRE_SPK); 
% discard spk ts at end of recording with inadequate post-spike period
spk_t = spk_t((spk_t+PST_SPK)*Fs < ecg_l); 

spk_n = length(spk_t);

% optional: use predetermined number of spikes
if spk_toggle
    if spk_n > NSPK
        spk_t(NSPK+1:end) = []; % only include the first NSPK spikes
        spk_n = NSPK;
    else
        error([nexfn ' does not exceed the # spike cutoff of' num2str(NSPK)...
            '.  Lower minimum cutoff for this unit.']);
    end
end

% create time vector
t = 1000*(-PRE_SPK:1/Fs:PST_SPK); % multiplied by 1000 to change to msec scale

% create time-reversed ecog for calculating confidence interval, using only
% the recording segment from which the spikes were taken
i_end = int32((spk_t(end)+PST_SPK)*Fs);
ecogallf=fliplr(ecogall(:,1:i_end));

% parse ecog centered at spike ts from each contact pair
E = zeros(spk_n,length(t),ecg_n); %3D matrix for storing spike triggered ecog data
Ef = zeros(spk_n,length(t),ecg_n); %3D matrix for storing spike triggered time-reversed ecog data

for i = 1:spk_n
    tmp1 = int32((spk_t(i)-PRE_SPK) * Fs);  % int32 used to keep index in integer format
    tmp2 = int32((spk_t(i)+PST_SPK) * Fs);
    % since int32 rounds to the next closest integer, there may be some
    % cases where the parsing indeces are off by 1.  fix them on the fly.
    d = tmp2 - tmp1;
    % lengt(t)-d must equal 1 for the parsed ecog data to fit E 3D matrix
    if length(t)-d == 0
        tmp2 = tmp2-1;
    elseif length(t)-d == 2
        tmp2 = tmp2+1;
    end
    for k = 1:ecg_n
        E(i,:,k) = ecogall(k,tmp1:tmp2);
        Ef(i,:,k) = ecogallf(k,tmp1:tmp2);
    end
end

STA = mean(E,1);  % spike-triggered ecog 
STA_mean = mean(STA,2); % mean of STA

STAf = mean(Ef,1);
% STAf_mean = mean(STAf,2); % mean of STAf
% STAf_std = std(STAf);

% % create randomly triggered ecog 
% tend = spk_t(end)-PST_SPK; % random timestamps are bound in time between 0 sec and time of last spike minus post-spike period
% rnd_t = ceil(1000*(tend*rand(spk_n,1))); % random time (msec)
% rnd_t = rnd_t/1000+PRE_SPK; %random time (sec), plus pre-spike period

% % parse ecog centered at spike ts from each contact pair
% Er = zeros(spk_n,length(t),ecg_n); %3D matrix for storing randomly triggered ecog data
% for i = 1:spk_n
%     tmp1 = int32((rnd_t(i)-PRE_SPK) * Fs);  % int32 used to keep index in integer format
%     tmp2 = int32((rnd_t(i)+PST_SPK) * Fs);
%     % since int32 rounds to the next closest integer, there may be some
%     % cases where the parsing indeces are off by 1.  fix them on the fly.
%     d = tmp2 - tmp1;
%     % lengt(t)-d must equal 1 for the parsed ecog data to fit E 3D matrix
%     if length(t)-d == 0
%         tmp2 = tmp2-1;
%     elseif length(t)-d == 2
%         tmp2 = tmp2+1;
%     end
%     for k = 1:ecg_n
%         Er(i,:,k) = ecogall(k,tmp1:tmp2);
%     end
% end

% RTA = mean(Er,1); % randomly triggered ecog average
% RTA_std = std(RTA); % standard deviation of RTA

%% plot figure
% hf1 = figure;
% hold on
% % figure 1: non-normalized, same plot
% 
% create a vector of ones used later for SD plotting
vones = ones(1,length(t));

% specify where along the time axis to place ecog pair label
t_text = -PRE_SPK*1e3*0.9;
% 
% % find min and max of STA
% vmin = min(min(min(STA)));
% vmax = max(max(max(STA)));
% 
% % calculate stacking constant used for plotting
% if (vmax-vmin) > 2*STDEV*max(STAf_std)
%     C_stk = vmax-vmin;
% else
%     %The 'else' case prevents stdev lines of neighboring ecog pairs from
%     %overlapping.
%     C_stk = 2*STDEV*max(STAf_std); 
% end
% 
% for k = 1:ecg_n
%     stk = (ecg_n - k) * C_stk;
%     z = STA(1,:,k) + stk; 
%     % create lines for STA average, upper and lower thresholds for plotting
%     avg = (STAf_mean(1,1,k)*vones) + stk;
%     thu = avg + (STDEV * STAf_std(1,1,k));
%     thl = avg - (STDEV * STAf_std(1,1,k));
%     plot(t,z,'LineWidth',2); % plot STA
%     plot(t,avg,'-k'); %plot time-reversed triggered average
%     plot(t,thu,'-.r',t,thl,'-.r'); % plot lower and upper thresholds
%     if k==M1_ch
%         text(t_text,thu(1),...
%             ['ecog' num2str(k) num2str(k+1) ',M1'],...
%             'VerticalAlignment','Top','FontWeight','bold'); % label M1 ecog pair
%     else
%         text(t_text,thu(1),...
%             ['ecog' num2str(k) num2str(k+1)],...
%             'VerticalAlignment','Top'); % label all other ecog pairs
%     end
% end
% xlabel('Time (millisec)');
% ylabel('Non-normalized STA (microV)');
% title([filename sprintf('\n')...
%     '# spks=' num2str(spk_n) ', sig thresh=' num2str(STDEV) 'SD'],...
%     'FontSize',12);
% ylm=ylim;
% plot([0 0],[ylm(1) ylm(2)],'k--');
% hold off


%% M1 STA processing
% 
% m1sta=STA(1,:,M1_ch);
% 
% % Build a finite impulse response (FIR)bandpass filter kernel with N "taps" 
% % and a user-defined bandpass setting.
% % (where sampling rate is 1000 Hz, nyquist frequency is
% %500 Hz, and the filter bandpass is 0.15 to 0.2 of the distance between 0 and 500).
% %the number of elements in the filter kernel is 100 to provide fast
% %roll-off.  The filter type is  windowed sinc filter, using the default
% %Hamming window
% 
% % Define the number of taps
% Ntap = 100;    
% % Define bandpass filter cutoffs, centered at peak M1 PSD
% [m1psd,freq]=pwelch(m1sta,WINDOW,NOVERLAP,NFFT,Fs);
% [max_pwr max_i]=max(m1psd);
% m1psdpeak=ceil(freq(max_i));
% COlo = m1psdpeak-2; % low cutoff of bandpass filter
% COhi = m1psdpeak+2; % high cutoff of bandpass filter
% 
% % Define a 2 element vector with the bandpass window
% WN = [COlo,COhi]/(Fs/2);
% % Build FIR bandpass filter
% B=fir1(Ntap,WN);
% % Plot the frequency response for bandpass filter kernel b to show that is
% % really is a bandpass filter 
% % freqz(B);
% 
% m1stafilt=conv(m1sta,B,'same');
% 
% % Hilbert transform
% xc = hilbert(m1stafilt);
% env=abs(xc);
% 
% % find max values for STA and 
% [vmax vi]=max(m1sta);
% vt=t(vi);
% [vfmax vfi]=max(env);
% vft=t(vfi);
% 
% % find STA at t=0
% vt0=m1stafilt(t==0);
% 
% % find phase
% % phase = atan(imag(xc)./real(xc));
% phase = atan2(imag(xc),real(xc));
% 
% % More plots
% hf2=figure;
% plot(t,m1sta),hold on
% plot(t,m1stafilt,'k','LineWidth',2)
% plot(t,env,'k');
% avg = STAf_mean(1,1,M1_ch)*vones;
% thu = avg + (STDEV * STAf_std(1,1,M1_ch));
% thl = avg - (STDEV * STAf_std(1,1,M1_ch));
% plot(t,thu,'-.r',t,thl,'-.r'); % plot lower and upper thresholds
% plot(vt,vmax,'o');
% plot(vft,vfmax,'ok');
% plot(0,vt0,'ok','MarkerFaceColor','k');
% ylm=ylim;
% plot([0 0],[ylm(1) ylm(2)],'k--'); % plot vertical bar for t=0
% xlabel('Time (millisec)');
% ylabel('Voltage (microV)');
% title([filename sprintf('\n')...
%     '# spks=' num2str(spk_n) ', sig thresh=' num2str(STDEV) 'SD'],...
%     'FontSize',12);
% text(t_text,thu(1),...
%     ['ecog' num2str(M1_ch) num2str(M1_ch+1) ', M1'],...
%     'VerticalAlignment','bottom',...
%     'FontWeight','bold'); % label all other ecog pairs
% text(250,thl(1),...
%     ['peak PSD @ ' num2str(freq(max_i)) 'Hz' sprintf('\n')...
%     'STA peak @ t=' num2str(vt) 'ms, max=' num2str(vmax) sprintf('\n')...
%     'env peak @ t=' num2str(vft) 'ms, max=' num2str(vfmax) sprintf('\n')...
%     'phase=' num2str(phase(t==0)) 'rad'],...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','top');
% legend('M1 STA',['M1 STA '  num2str(COlo) '-' num2str(COhi) 'Hz filtered'],'amplitude envelope');
% legend('boxoff');
% 
% hold off

% More plots
hf3=figure;

m1sta=STA(1,:,M1_ch);
m1staf=STAf(1,:,M1_ch);

STAf_mean = mean(STAf(:,401:600,:),2); % mean of STAf
STAf_std = std(STAf(:,401:600,:));

plot(t,m1sta,'k','LineWidth',2),hold on
plot(t,m1staf);

avg=STAf_mean(1,1,M1_ch)*vones;

thu2 = avg + (2 * STAf_std(1,1,M1_ch));
thl2 = avg - (2 * STAf_std(1,1,M1_ch));
thu225 = avg + (2.25 * STAf_std(1,1,M1_ch));
thl225 = avg - (2.25 * STAf_std(1,1,M1_ch));
thu25 = avg + (2.5 * STAf_std(1,1,M1_ch));
thl25 = avg - (2.5 * STAf_std(1,1,M1_ch));
thu275 = avg + (2.75 * STAf_std(1,1,M1_ch));
thl275 = avg - (2.75 * STAf_std(1,1,M1_ch));
thu3 = avg + (3 * STAf_std(1,1,M1_ch));
thl3 = avg - (3 * STAf_std(1,1,M1_ch));

plot(t,avg,'-k'); %plot time-reversed triggered average

plot(t,thu2,'-.r',t,thl2,'-.r'); % plot lower and upper thresholds
plot(t,thu225,'-.r',t,thl225,'-.r'); % plot lower and upper thresholds
plot(t,thu25,'-.r',t,thl25,'-.r'); % plot lower and upper thresholds
plot(t,thu275,'-.r',t,thl275,'-.r'); % plot lower and upper thresholds
plot(t,thu3,'-.r',t,thl3,'-.r'); % plot lower and upper thresholds

ylm=ylim;
plot([0 0],[ylm(1) ylm(2)],'k--'); % plot vertical bar for t=0
xlabel('Time (millisec)');
ylabel('Voltage (microV)');
title([filename sprintf('\n')...
    '# spks=' num2str(spk_n) ', SD calculated between -100 and +100 msec'],...
    'FontSize',12);
text(t_text,thu2(1),...
    ['ecog' num2str(M1_ch) num2str(M1_ch+1) ', M1'],...
    'VerticalAlignment','bottom',...
    'FontWeight','bold'); % label all other ecog pairs
text(400,thu2(1),'2SD',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(400,thu225(1),'2.25SD',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(400,thu25(1),'2.5SD',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(400,thu275(1),'2.75SD',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
text(400,thu3(1),'3SD',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
legend('spike-triggered','time-reversed');
legend('boxoff');

%% save variables and figures
% save([filename '_sta'],'t','STA','STA_mean',...
%     'm1sta','m1stafilt','xc','env','vmax','vt','vfmax','vft','vt0','phase');
% saveas(hf1,[filename '_sta'],'fig');
% saveas(hf2,[filename '_m1sta'],'fig');


% figure 2: normalized, same plot

% hold on
% for k = 1:ecg_n
%     vmin = min(STA(1,:,k)); % max value of average
%     vmax = max(STA(1,:,k)); % min value of average
%     C = ecg_n - k; % constant added to stack the waves for comparison such that the first contact pair is plotted at the top
%     z = (STA(1,:,k)-vmin)/(vmax-vmin) + C; 
%     plot(t,z);
% end
% xlabel('Time (msec)');
% ylabel('normalized ecog');
% title(filename);
% ylm = ylim;
% plot([0 0],[ylm(1) ylm(2)],'k--')
% hold off

% figure 3: non-normalized, subplots

% for k=1:ecg_n
%     subplot(ecg_n,1,k)
%     plot(t,STA(1,:,k));
%     ylim([-3 2]);
%     if k==1
%         title(filename);
%         xlabel('Time (msec)');
%         ylabel('ecog');
%     end
% end



