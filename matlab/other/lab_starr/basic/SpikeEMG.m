function SpikeEMG()
% SpikeEMG()
% This program analyses Spike times versus EMG data
% Created by RST, 2002-04-05

global RTPATH	% starting directory containing NEX file
global fs       % Sampling rate for processed EMG & spike SDF
fudge = 5;      % Conversion to make sdf = sp/sec
fs = 1000;      % sampling rate for EMG & SDF
p = path;

[nexname, RTPATH] = uigetfile('*.nex', 'Select spike time file (NEX)');

if (nexname ~= 0)
    cd(RTPATH);
    [nvar, varname, types] = nex_info(nexname);
    if nvar > 1
        error([num2str(nvar) ' spikes in ' nexname '!!!  I don''t know how to process > 1 spike yet']);
    end
    [spk.n, spk.t] = nex_ts(nexname,varname);
else
    error(['I can''t find the NEX file:  ' nexname ' in ' RTPATH]);
end

fname = strrep(nexname,'.nex','');

% Output a txt file for Labview or Matlab ISI analyses
write_isi_txt(spk.t,fname);

tmp = dir([fname '.mat']);
[sz x] = size(tmp);
if(sz ==0)
    error(['I can''t find the EMG file:  ' tmp ' in ' RTPATH]);
end
load( [fname '.mat'] );
[emg_n, data_len] = size(emg_chan);

% Convert spike times to delta function
spk.delt = spk_t2delta(spk.t,data_len);

% convolve with gaussian to make spike density function
load gaus_flt gaus_20ms_1kHz
flt_len = (length(gaus_20ms_1kHz)-1)/2;
spk.sdf = fudge*conv(gaus_20ms_1kHz,spk.delt);
spk.sdf(1:flt_len)=[];
spk.sdf(data_len+1:data_len+flt_len)=[];

% Set up to plot data
left = 0.05;
width = .5-left;
height = 0.90/(n_emg + 1);      % Used for both time & coherence plots
figure

% Plot time domain
bottom = 0.03-height;
for i = emg_n:-1:1
    bottom = bottom + height;
    subplot('position',[left bottom width height]);
    plot(time,emg_chan(i,:),'k-');
    axis([ min(time) max(time) 0 1.1*max(emg_chan(i,:))]);
    ylabel([ 'EMG-', num2str(i)]);
    if i == emg_n
        set(gca,'ytick',[ 0 ]);
        xlabel('Time (sec)');
    else
        set(gca,'ytick',[ 0 ],'xtick',[ ]);
    end
end
bottom = 1-height;
subplot('position',[left bottom width height]);
plot(time,spk.sdf,'k-');
axis([ min(time) max(time) 0 1.1*max(spk.sdf)]);
xlabel('Time (sec)');
ylabel('Spikes/sec');
hold on
% Size to make it look good
c = get(gcf);
c.Position(2) = 275;
c.Position(3) = 870;
c.Position(4) = 680;
set(gcf,'Position',c.Position);


%-------------------------------------
% For testing !!! Introduce sine wave
%-------------------------------------
freq = 8;
sm_sin = 10*sin(freq*time*2*pi);
sm_sin(2,:) = sm_sin;
sm_sin(3,:) = sm_sin(1,:);
sm_sin(4,:) = sm_sin(1,:);
%emg_chan = emg_chan + sm_sin;
%spk.sdf = spk.sdf + sm_sin(1,:);

[spk.auto, spk.autot] = autocorr(1000*diff(spk.t),500,10);
%-------------------------------------
% Compute SDF-EMG cross-correlations
%-------------------------------------


nfft = 1024;        % Size of windows
overlap = nfft-round(0.1*nfft); % Overlap windows
cutoff = 50;        % Show only frequencies < cutoff
coh = zeros(emg_n, (nfft/2+1) );

bottom = 0.03-height;
left = 0.8;
width = 0.98-left;
for i = emg_n:-1:1
    bottom = bottom + height;
    subplot('position',[left bottom width height]);

    [coh,f] = cohere(squeeze(spk.sdf),squeeze(emg_chan(i,:)),...
        nfft,fs,[],overlap,'linear');
    sel = find(f<cutoff);
    plot(f(sel),coh(sel),'k-');
    if i == emg_n
        xlabel('Coherence (Hz)');
    end

    % Title plot
    text(cutoff-13,0.8,['\bfPow-', num2str(i)]);
    set(gca,'ytick',[ 0.2 0.4 0.6 0.8 ],'xtick',[ 20 40 ]);
    axis( [ 0 cutoff 0 1 ] );
end

