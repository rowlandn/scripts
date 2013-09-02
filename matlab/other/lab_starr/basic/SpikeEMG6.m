function [spk,f_sum,cl] = SpikeEMG6(nexname)
% [spk,emg] = SpikeEMG6(fname)
% This program analyses Spike times versus EMG data
% Created by RST, 2002-04-05
% Modified by RST, 2003-10-15
% Modified by RST, 2004-01-27
%	

global RTPATH	% starting directory containing NEX file
global FS       % Sampling rate for processed EMG & spike SDF
global FSS		% Subsample SDF & EMG to this rate
global SRCH_LO	% low freq for search for significant oscillations
global SRCH_HI	% top freq for search for significant oscillations
global SIG_SNR	% significant oscill must have power SNR > SIG_SNR * std(power)
global SIG_OSC	% or OSC_IND > SIG_OSC
global SIG_COH
global MAX_N	% find max's in acorr over MAX_N points
global AUTOCORR_LAG
global SEG_PWR

NFFT = 2^SEG_PWR;	% For spectral analysis of spike autocorr
SS_INTERV = 1000/FSS;

% Constants for Plotting xcorr & coh
MX_FRQ = 30;	% for plotting

global FRQ_COL	% coherence in 4th column
global COH_COL	% coherence in 4th column
global PHS_COL	% phase in 5th column

TOP = 0.97;		% Top margin of page
LEFT = 0.05;	% Left margin of page

% Load filter information
load('abfconv_filters');

if (nexname ~= 0)
    [nvar, varname, types] = nex_info(nexname);
    if nvar > 1
        error([num2str(nvar) ' spikes in ' nexname '!!!  I don''t know how to process > 1 spike yet']);
    end
    [spk.n, spk.t] = nex_ts(nexname,varname);
else
    error(['I can''t find the NEX file:  ' nexname ]);
end

fname = strrep(nexname,'.nex','');

tmp = dir([fname '.mat']);
[sz x] = size(tmp);
if(sz ==0)
    warning(['I can''t find the EMG file:  ' fname ]);
	EMG_FLAG = false;
	n_emg = 1;
else
	EMG_FLAG = true;
	load( [fname '.mat'] );		% NOTE:  emgs come in @ 1kHz samped w/ LP filter @ 100 Hz
	[n_emg, data_len] = size(emg_chan);
end

% Convert spike times to delta function & sdf
if EMG_FLAG
	spk.delt = spk_t2delta(spk.t,data_len);		
else
	spk.delt = spk_t2delta(spk.t);
	data_len = length(spk.delt);
	time = (1/FS) : 1/FS : (data_len/FS);
end
temp = spk2sdf(spk.delt);		%NOTE: spk2sdf convolves w/ 10 ms gaussian, little power above 100 Hz
spk.sdf = filtfilt(LowPass_100Hz_1kHz.tf.num, 1, temp ); % Chop out > 100 Hz 
spk.sdf = decim_rst(spk.sdf,SS_INTERV);		% subsample to 200 Hz
time5 = decim_rst(time,SS_INTERV);
spk.dur = max(time);

% Compute autocorr & spectrum from spike times
[spk.acorr,spk.lag] = Spike_autocorr(spk.t,AUTOCORR_LAG);
% find & remove central valley
avg_autocorr = mean(spk.acorr(1:(AUTOCORR_LAG-10)));
for i = 1:AUTOCORR_LAG
	if spk.acorr(i+AUTOCORR_LAG+1) >= avg_autocorr
		break
	end
end
spk_corr = spk.acorr;
for j = -i:i
	spk_corr(j+AUTOCORR_LAG+1) = avg_autocorr;
end
spk_corr = detrend(spk_corr);
temp = filtfilt(LowPass_100Hz_1kHz.tf.num, 1, spk_corr ); % Chop out > 100 Hz 
spk.acorr_subsamp = decim_rst(temp,SS_INTERV);	% sub-sample autocorr down to 200 Hz
spk.lag_subsamp = decim_rst(spk.lag,SS_INTERV);	% sub-sample lags down to 200 Hz

% Do spectral computation @ 200 Hz 
[spk.pow,spk.frq] = periodogram_power(spk.acorr_subsamp,[],NFFT,FSS);
spk.acorr(1:AUTOCORR_LAG/SS_INTERV) = [];
spk.lag(1:AUTOCORR_LAG/SS_INTERV) = [];

% Find peaks in power spectrum 
srch_ind = find( spk.frq > SRCH_LO & spk.frq < SRCH_HI );
pow_pks = find_local_nmax(spk.pow, MAX_N);
[spk.pk_ind, ind_pow_pks, ind_srch_ind] = intersect( pow_pks, srch_ind );

% Compute SNR for each peak
spk.snr = (spk.pow(spk.pk_ind) - mean( spk.pow(srch_ind))) ./ std(spk.pow) ;

% Compute Oscil index for each peak
pow_area = find_local_areas( spk.pow, mean(spk.pow(srch_ind)), pow_pks );
spk.oscil = 100 .* pow_area(ind_pow_pks)' ./ sum(spk.pow);		% Oscillation index

% Find significant peaks for plotting
a = find(spk.snr > SIG_SNR);
b = find(spk.oscil>SIG_OSC);
spk.sig_pk_ind = union( a, b);
spk.sig_pks = spk.pk_ind( spk.sig_pk_ind );

%%%%%%%%%%%%%% Unit plotting
% Set up to plot time data
figure
set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
left = LEFT;
width = 0.45-left;
height = 0.90/(n_emg + 1);      % Used for both time & coherence plots

bottom = TOP-height;
subplot('position',[left bottom width height]);
plot(time5,spk.sdf,'k-');
axis([ min(time) max(time) 0 1.1*max(spk.sdf)]);
xlabel('Time (sec)');
ylabel('Spikes/sec');
title([nexname]);
set(gca,'FontSize',8);
% Size to make it look good
c = get(gcf);
c.Position(2) = 275;
c.Position(3) = 870;
c.Position(4) = 680;
set(gcf,'Position',c.Position);

% Setup to plot spectrum & autocorr
left = left+width+0.05;
width = 0.95-left;
height = 0.90/(n_emg + 1);      % Used for both time & coherence plots
bottom = TOP-height;

% Plot spectra & autocorr
subplot('position',[left bottom width height]);
plot(spk.frq,spk.pow,'k-','LineWidth',2);
hold on
plot(spk.frq(spk.sig_pks),spk.pow(spk.sig_pks),'ro','LineWidth',2);
xlim([0 MX_FRQ]);
ax1 = gca;
set(ax1,'Box','off','FontSize',8);
xtks = get(ax1,'XTick');
xtks(:,end) = [];
set(ax1,'XTick',xtks);
ytks = get(ax1,'YTick');
ytks(:,end) = [];
set(ax1,'YTick',ytks);
ylabel('Power');

ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','XColor','b',...
	'YAxisLocation','right','Color','none','YColor','b');
	
% Plot autocorr on second axis
h2 = line(spk.lag_subsamp,spk.acorr_subsamp,...
	'Color','b','Parent',ax2);
axis([0,max(spk.lag_subsamp),min(spk.acorr_subsamp),max(spk.acorr_subsamp)]);

xtks = get(ax2,'XTick');
xtks(:,1) = [];
xtks(:,end) = [];
set(ax2,'XTick',xtks,'FontSize',8);
ytks = get(ax2,'YTick');
ytks(:,1) = [];
set(ax2,'YTick',ytks);
xlabel('Lag (ms)');
%%%%%%%%%%%%%% END plotting of unit data

if EMG_FLAG
	% Sub-sample EMG data
	emg = zeros( n_emg,length(time5) );
	for i = 1:n_emg
		emg(i,:) = decim_rst(emg_chan(i,:),SS_INTERV);
	end

	% Setup to Compute xcorr & coherence between spike times & emg
	[f,t,cl] = sp2a2_m(spk.sdf',emg(1,:)',FSS,SEG_PWR,'t');
	f_sum = zeros(length(f),5,n_emg);
	t_sum = zeros(length(t),2,n_emg);
	cl.coh_sig = 1 - SIG_COH^(1/(cl.seg_tot-1));	% From Halliday
	cl.what = sprintf('Spike-vs-EMG/ACCEL');	% Label for coh subplots
	freq_pts=round(MX_FRQ/cl.df);
	ch_max = 0;
	
	% Compute xcorr & coherence between spikes & emg
	for i = 1:n_emg
		[f,t,tmp] = sp2a2_m(spk.sdf',emg(i,:)',FSS,SEG_PWR,'t');
		mx = max( f(1:freq_pts,COH_COL) );	% Find max coherence across estimates
		if mx > ch_max		ch_max = mx;	end
		% Convert phase to time-lag/lead
		f(:,PHS_COL) = -1000 .* f(:,PHS_COL)./((2*pi) .* f(:,FRQ_COL)) ;
		f_sum(:,:,i) = f;
		t_sum(:,:,i) = t;
	end
	
	%%%%%%%%%%%%%% Plotting EMG data
	% Set up to plot time data
	left = LEFT;
	width = 0.45-left;
	height = 0.90/(n_emg + 1);      % Used for both time & coherence plots
	bottom = 0.03-height;
	% Plot EMG time domains
	for i = n_emg:-1:1
        bottom = bottom + height;
        subplot('position',[left bottom width height]);
        plot(time5,emg(i,:),'k-');
        if i < n_emg 
			ylabel([ 'EMG-', num2str(i)]);
            set(gca,'ytick',[ 0 ],'xtick',[ ]);
			axis([ min(time) max(time) 0 1.1*max(emg(i,:))]);
        elseif i == n_emg & i<5
			ylabel([ 'EMG-', num2str(i)]);
			axis([ min(time) max(time) 0 1.1*max(emg(i,:))]);
            set(gca,'ytick',[ 0 ]);
            xlabel('Time (sec)');
        elseif i == n_emg & i==5
			ylabel('ACCEL');
			axis([ min(time) max(time) 1.1*min(emg(i,:)) 1.1*max(emg(i,:))]);
            set(gca,'ytick',[ 0 ]);
            xlabel('Time (sec)');
        end
		set(gca,'FontSize',8);
	end
	% Setup to plot coherences
	left = left+width+0.05;
	width = 0.95-left;
	height = 0.90/(n_emg + 1);      % Used for both time & coherence plots
	
	if ch_max < (2.*cl.coh_sig)	ch_max = (2.*cl.coh_sig);	end
	
	% Plot spike/EMG coherences
	bottom = 0.03-height;
	for i = n_emg:-1:1
        bottom = bottom + height;
        subplot('position',[left bottom width height]);
		psp_coh_phase(f_sum(:,:,i),cl,MX_FRQ,ch_max);
	end
	%%%%%%%%%%%%%% END of Plotting EMG data
else
	f_sum = NaN;
	cl = NaN;
end
