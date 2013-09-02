function spk = SpikeAutocorr(nexname, varnum)
% spk = SpikeAutocorr(fname, varnum)
% This program analyses Spike autocorrelations from NEX files
% Created by RST, 2002-04-05
% Modified by RST, 2003-10-15
% Modified by RST, 2004-01-27
% Modified by RST, 2005-07-12 to pick indicated variable number
%
% Optional argument 'varnum' can be used to process matching variable
% number

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

% Constant for Plotting
TOP = 0.97;		% Top margin of page
LEFT = 0.05;	% Left margin of page
MX_FRQ = 30;	% for plotting

% Variables used in Rory's Nexload, added here for compatibility
SPIKE_TYPE = 0;
n = 1;
isi_batching = 0;	% used for something in Nexload, added here merely for compatibility
Nexdir(n).name = nexname;	% More compatibility w/ Nexload
nexpath = '.\';
verbose = 0;

% Load smoothing filter
load('abfconv_filters');

if exist(nexname) ~= 2
    error(['I can''t find the NEX file:  ' nexname ]);
else
	Nexload		% Run Rory's script for reading a nex file
end

spk.name = deblank( strrep(nexname,'.nex',['_' Nexvar.Name(varnum,:)] ) );
spk.t = Nexdata.TS{ varnum } ;
spk.n = Nexdata.TS_count( varnum );

% Convert spike times to delta function & sdf
spk.delt = spk_t2delta(spk.t);
data_len = length(spk.delt);
time = (1/FS) : 1/FS : (data_len/FS);
temp = spk2sdf(spk.delt);		%NOTE: spk2sdf convolves w/ 10 ms gaussian, little power above 100 Hz
spk.sdf = filtfilt(LowPass_100Hz_1kHz.tf.num, 1, temp ); % Chop out > 100 Hz 
spk.sdf = decim(spk.sdf,SS_INTERV);		% subsample to 200 Hz
time5 = decim(time,SS_INTERV);
spk.dur = max(time);

% Compute autocorr & spectrum from spike times
[spk.acorr,spk.lag] = Spike_autocorr(spk.t,AUTOCORR_LAG);

% Find significant oscillations in autocorr
[spk.snr, spk.oscil, spk.pow, spk.frq, spk.pk_ind, ...
	spk.acorr_subsamp, spk.lag_subsamp] = autocorr_anal(spk.acorr, spk.lag, FSS);
	
% Find significant spectral peaks for plotting
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
