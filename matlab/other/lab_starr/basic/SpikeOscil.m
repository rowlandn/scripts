function spk = SpikeOscil( spk_t, unitname, show)
% spk = SpikeOscil(spk.t, unitname)
% This program searches for significant oscillations in spiketrains
% Inputs:
%	spk.t = vector times of spikes (in seconds)
%	unitname - string containing cat'd fname & unitname
%	show - optional flag to control display of graphical output
%		(default:true)
% Outputs:
%	spk = structure containing results
%
% Created by RST, 2005-07-30
%

global SEG_PWR
global FS       % Sampling rate
global SIG_OSC	% threshold for significance
global CNTL_FRQ	% Frq range used to compute spectral SD
global N_SHUF	% number of ISI shuffles for control spectra
global AUTOCORR_LAG
global SRCH_LO	% low freq for search for significant oscillations
global SRCH_HI	% top freq for search for significant oscillations
global MAX_N	% find max's in spectra over MAX_N symmetric points (must be odd) 7 => 3.5 Hz separation

if ~exist('show','var')
	show = true;
end

% Load filter information
load('abfconv_filters');

NFFT = 2^SEG_PWR;	% For spectral analysis of spike autocorr
spk.name = unitname;

% Compute autocorr 
[spk.acorr,spk.lag] = Spike_autocorr(spk_t,AUTOCORR_LAG);
spk.acorr(1:(AUTOCORR_LAG-1)) = [];
spk.lag(1:(AUTOCORR_LAG-1)) = [];
acorr_lp = filtfilt(LowPass_100Hz_1kHz.tf.num, 1, spk.acorr ); % Chop out > 100 Hz 


% Convert spike times to delta function
spk.t = round( 1000.*spk_t );
len_delt = spk.t(end)-spk.t(1)+1;
spk.delt = zeros(1,len_delt);
spk.delt( spk.t - spk.t(1)+1 ) = 1;

% Calculate the power on the spike train
hann = hanning(NFFT);
[ pow, freq] = psd( spk.delt-mean(spk.delt), NFFT, FS, hann); 
spk.pow = pow';
spk.freq = freq';
len_spect = length(freq);

% Calculate psd's from globaly shuffled spike trains with the same ISI
fprintf('Shuffling');
clk_tk = round(N_SHUF/10);
isi = diff(spk.t);
isi = isi(isi>0); %make sure isi's are good
pow_shuf=zeros( N_SHUF, len_spect);
for i=1:N_SHUF
    r=randperm(length(isi));
    rand_isi(r)=isi;
    y=cumsum(rand_isi);
    rand_delt = zeros(1,len_delt);
    rand_delt([1 y+1])=1; 
    [pow_rand,freq_rand] = psd( rand_delt-mean(rand_delt), NFFT, FS, hann);
    pow_shuf(i,:)=pow_rand';
	if ~mod(i,clk_tk);	fprintf('.');	end
end
spk.pow_rand = mean(pow_shuf);	% Mean pow of randomized ISI's
fprintf('\n');

% Now compute ISI-compensated normalized spectrum
spk.pow_comp = spk.pow ./ spk.pow_rand;
% spk.pow_comp = spk.pow; %3/2/09:uncomment this line and comment out line above to use non-normalized spectra

% Compute threshold for significance
cntl_ind = find(freq>=CNTL_FRQ(1) & freq<CNTL_FRQ(2) );
p = SIG_OSC/len_spect;
spk.sig_thresh = norminv(1-p) * std( spk.pow_comp(cntl_ind) ) + 1;
% spk.sig_thresh = norminv(1-p) * std( spk.pow_comp(cntl_ind) ) + mean( spk.pow_comp(cntl_ind) ); %3/2/09: uncomment this line for non-normalized spectra

% find significant peaks w/in search range
srch_inds = find( spk.freq>=SRCH_LO & spk.freq<SRCH_HI );
spk.sig_inds = find_sig_peaks( spk.pow_comp, ...
	spk.sig_thresh, MAX_N, srch_inds );

if show
	%%%%%%%%%%%%%% Plotting
	% Set up axes
	MARGIN = 0.06;	
	TOP = 1-MARGIN;		% Top margin of page
	LEFT = MARGIN;	% Left margin of page
	ZOOM = 0.15;		% Fraction of plots to magnify in Left column of plots
	ZOOM_WID = 0.3;
	HEIGHT = 0.4;	% Height of plots, used for both spect & autocorr
	figure
	set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
	% Size to make it look good
	c = get(gcf);
	c.Position(2) = 275;
	c.Position(3) = 870;
	c.Position(4) = 680;
	set(gcf,'Position',c.Position);

    % Plot Zoomed oscil results
    left = MARGIN;
    width = ZOOM_WID;
    height = HEIGHT;      % Used for both spectrum & autocorr plots
    bottom = TOP-height;
    subplot('position',[left bottom width height]);
    plot(spk.freq, spk.pow_comp,'k-','LineWidth',1);
    xlim([0 SRCH_HI*ZOOM]);
    hold on
    plot(spk.freq, spk.sig_thresh,'r:','LineWidth',2);
    sig_inds = spk.sig_inds;
    plot(spk.freq(sig_inds), spk.pow_comp(sig_inds),'ro','LineWidth',2);
    xlabel('Frequency (Hz)');
    ylabel('Normalized power');
    set(gca,'FontSize',8);

	% Plot Zoomed autocorr results
	bottom = bottom-height-MARGIN;
	subplot('position',[left bottom width height]);
	plot(spk.lag, spk.acorr,'k-','LineWidth',1);
	xlim([0 50]);
	lim = ylim;
	ylim([0 lim(2)]);
	xlabel('Lag (msec)');
	ylabel('spikes/sec');
	set(gca,'FontSize',8);

	% Plot whole oscil results
	left = ZOOM_WID+2*MARGIN;
	width = 1-(left+MARGIN);
	height = HEIGHT;      % Used for both spectrum & autocorr plots
	bottom = TOP-height;
	subplot('position',[left bottom width height]);
	plot(spk.freq, spk.pow_comp,'k-','LineWidth',1);
    xlim([0 SRCH_HI]);
	hold on
	plot(spk.freq, spk.sig_thresh,'r:','LineWidth',2);
	sig_inds = spk.sig_inds;
	plot(spk.freq(sig_inds), spk.pow_comp(sig_inds),'ro','LineWidth',2);
	xlabel('Frequency (Hz)');
	set(gca,'FontSize',8);
	title([ unitname ],'Interpreter','none','FontSize',14);

	% Plot whole autocorr results
	bottom = bottom-height-MARGIN;
	subplot('position',[left bottom width height]);
	plot(spk.lag, acorr_lp,'k-','LineWidth',1);
	xlim([0 spk.lag(end)]);
	ymin = min(acorr_lp(10:end))-10;
	ymax = max(acorr_lp(10:end));
	ylim([ymin ymax]);
	xlabel('Lag (msec)');
	set(gca,'FontSize',8);
end
