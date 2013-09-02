function spk = SpikeAutoAlex(ts, fname)
% [spk] = SpikeAutoAlex(timestamp)
% This program analyses Spike time autocorrelation
% Input time stamps in sec.s
% Created by RST, 2002-04-05
% Modified by RST, 2003-10-15
% Modified by RST, 2004-01-27
%	
SEG_PWR = 9;	% Segment length - specified as power of 2
FS = 1000;
FSS = 200;		% Subsample to this rate
SRCH_LO = 2.0;	% low freq for search for significant oscillations in autocorr
SRCH_HI = 30;	% top freq for search for significant oscillations
SIG_SNR = 7;	% significant oscill must have power SNR > SIG_SNR * std(power)
SIG_OSC = 10;	% or OSC_IND > SIG_OSC
SIG_COH = 0.0005;	% Significance threshold for Spike-EMG coherence :0.0005 => false positive in 5% of cases

MAX_N = 9;		% find max's in acorr over MAX_N symmetric points (must be odd) 9 => 3.5 Hz separation
AUTOCORR_LAG = 500;

% Indices into columns of Spike/EMG data (returned from NeuroSpec)
FRQ_COL = 1;	% frequency in 1st column
COH_COL = 4;	% coherence in 4th column
PHS_COL = 5;	% phase in 5th column

N_ACorr = 2;	% max # of significant autocorr pks reported in txt output
N_XCorr = 5;	% max # of significant spike/EMG coh pks reported in txt output

COH_FBAND = [ 0 1 ];	% Do XCorr analysis for sig cohs in this band
XCORR_LAG = 4000;		% MSec lags/lead for xcorrelation

NFFT = 2^SEG_PWR;	% For spectral analysis of spike autocorr
SS_INTERV = 1000/FSS;

% Constants for Plotting xcorr & coh
MX_FRQ = 30;	% for plotting

TOP = 0.97;		% Top margin of page
LEFT = 0.05;	% Left margin of page

% Load filter information
load('abfconv_filters');
spk.t = ts;
	spk.delt = spk_t2delta(spk.t);
	data_len = length(spk.delt);
	time = (1/FS) : 1/FS : (data_len/FS);

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
height = 0.90;      % Used for both time & coherence plots

bottom = TOP-height;
subplot('position',[left bottom width height]);
plot(time5,spk.sdf,'k-');
axis([ min(time) max(time) 0 1.1*max(spk.sdf)]);
xlabel('Time (sec)');
ylabel('Spikes/sec');
if exist('fname')
    title(fname);
end
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
height = 0.90;      % Used for both time & coherence plots
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

