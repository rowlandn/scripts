function [spk,f_sum,cl_sum] = SpikeEMG5(fname)
% [spk,emg] = SpikeEMG5(fname)
% This program analyses Spike times versus EMG data
% Created by RST, 2002-04-05
% Modified by RST, 2003-10-15
% Modified by RST, 2004-01-27
%	

global FS       % Sampling rate for processed EMG & spike SDF
global SRCH_LO
global SRCH_HI
global SIG_SNR
global SIG_OSC
global SIG_COH
global FRQ_COL
global COH_COL
global PHS_COL
global MAX_N
global AUTOCORR_LAG
global SEG_PWR
global PLOTTING

NFFT = 2^SEG_PWR;	% For spectral analysis of spike autocorr

% Constants for Plotting xcorr & coh
MX_FRQ = 30;	% for plotting

TOP = 0.97;		% Top margin of page
LEFT = 0.05;	% Left margin of page

% Load filter information
load('abfconv_filters');

if isempty(fname)
	[fname, pathname, filterindex] = uigetfile('*.mat','Select test file to process');
end

tmp = dir(fname);
[sz x] = size(tmp);
if(sz ==0)
    error(['I can''t find the file:  ' fname ]);
else
	EMG_FLAG = true;
	load( fname );		% NOTE:  emgs come in @ 1kHz samped w/ LP filter @ 100 Hz
	[n_emg, data_len] = size(emg_chan);
end

temp = spk2sdf(spk.delt);		%NOTE: spk2sdf convolves w/ 10 ms gaussian, little power above 100 Hz
spk.sdf = filtfilt(LowPass_100Hz_1kHz.tf.num, 1, temp ); % Chop out > 100 Hz 
spk.sdf = decim(spk.sdf,5);		% subsample to 200 Hz
time5 = decim(time,5);
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
spk.acorr_subsamp = decim(temp,5);	% sub-sample autocorr down to 200 Hz
spk.lag_subsamp = decim(spk.lag,5);	% sub-sample lags down to 200 Hz

% Do spectral computation @ 200 Hz 
[spk.pow,spk.frq] = periodogram_power(spk.acorr_subsamp,[],NFFT,200);
spk.acorr(1:AUTOCORR_LAG/5) = [];
spk.lag(1:AUTOCORR_LAG/5) = [];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process EMG data
% Sub-sample EMG data
emg = zeros( n_emg,length(time5) );
for i = 1:n_emg
	emg(i,:) = decim(emg_chan(i,:),5);
end

% Setup to Compute xcorr & coherence between spike times & emg
[f,t,cl] = sp2a2_m(spk.sdf',emg(1,:)',200,SEG_PWR,'t');
f_sum = zeros(length(f),5,n_emg);
t_sum = zeros(length(t),2,n_emg);
cl_sum = cl;
cl_sum.coh_sig = 1 - SIG_COH^(1/(cl.seg_tot-1));	% From Halliday
freq_pts=round(MX_FRQ/cl.df);
ch_max = 0;

% Compute xcorr & coherence between spikes & emg
for i = 1:n_emg
	[f,t,cl] = sp2a2_m(spk.sdf',emg(i,:)',200,SEG_PWR,'t');
	mx = max( f(1:freq_pts,COH_COL) );	% Find max coherence across estimates
	if mx > ch_max		ch_max = mx;	end
	% Convert phase to time-lag/lead
	f(:,PHS_COL) = 1000 .* f(:,PHS_COL)./((2*pi) .* f(:,FRQ_COL)) ;
	f_sum(:,:,i) = f;
	t_sum(:,:,i) = t;
end
		

if PLOTTING
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
		title(fname);
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
		
		if ch_max < (2.*cl_sum.coh_sig)	ch_max = (2.*cl_sum.coh_sig);	end
		
		% Plot spike/EMG coherences
		bottom = 0.03-height;
		for i = n_emg:-1:1
            bottom = bottom + height;
            subplot('position',[left bottom width height]);
			psp_coh_phase(f_sum(:,:,i),cl_sum,MX_FRQ,ch_max);
		end
	%%%%%%%%%%%%%% END of Plotting EMG data
end
