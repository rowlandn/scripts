function [pk_val, pk_lag, xcor_sum] = SpikeEMGXCorr2(spk,emg,chan_sub,xcor_coh,xcor_frq)
%  [pk_val, pk_lag] = SpikeEMGXCorr(spk,emg,chan_sub,xcor_frq)
% This program analyses Spike times versus EMG data
% Created by RST, 2002-04-05
% Modified by RST, 2003-10-15
% Modified by RST, 2004-01-27
% Modified by RST, 2007-02-23  -  spike & emg data received rather than
% re-read from files
%	
%	fname = name of nex file - must be accompanied by matching mat file
%	chan_sub = EMG channel numbers w/ significant coherence in band of
%		interest
%	xcor_frg = frequencies of significant coherences
%
%	Returns:	pk_val = vector of peak R values for each sig coherence
%				pk_lag = time lag/lead for peak R value
%				xcor_sum = total of cross correlation
%

global RTPATH	% starting directory containing NEX file
global FS       % Sampling rate for processed EMG & spike SDF
global FSS       % Sampling rate for processed EMG & spike SDF
global SRCH_LO	% Used to set low frq range for band search
global SRCH_HI	% Used to set high frq range for band search
global SIG_COH
global MAX_N	% find max's in acorr over MAX_N points
global XCORR_LAG
global COH_FBAND

SS_INTERV = 1000/FSS;
FILT_LEN = 500;	% For bandpass filter
% Constants for Plotting xcorr & coh
MX_FRQ = 30;	% for plotting

global FRQ_COL	% coherence in 4th column
global COH_COL	% coherence in 4th column
global PHS_COL	% phase in 5th column

TOP = 0.97;		% Top margin of page
LEFT = 0.06;	% Left margin of page
TXT_MRG = 0.02;	% Margin for text

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process significant cohs in range, perform bandpass
	chan_sub = unique(chan_sub);
	n_xcorr = length(chan_sub);
	% Make bandpass filter
	xcorr_lag = round(XCORR_LAG * FSS / 1000);		% Lags in samples
	xcorr_tlag = -XCORR_LAG:SS_INTERV:XCORR_LAG;	% Time lags for plotting
	data_len = length(spk.sdf);
	[d,n_emg] = size(emg);
	if COH_FBAND(1) == 0;
		Wn = COH_FBAND(2);
	else
		Wn = COH_FBAND;
	end
	Wn = Wn .* 2 ./ FSS;	% Express cutoffs in terms of Nyquist limit
	b = fir1(FILT_LEN,Wn);	% Make filter
	
	% Bandpass Filter sdf & emgs
	bp.sdf = filtfilt_rst(b, 1, spk.sdf ); % Band-pass SDF to frq range of interest
	bp.sdf = detrend(bp.sdf);
	bp.sdf = bp.sdf ./ (max(bp.sdf)-min(bp.sdf));	% Normalize spike data

	bp.emg = zeros(length(chan_sub),data_len);
	bp.xcor = zeros(length(chan_sub),2*xcorr_lag+1);
	for i=1:length(chan_sub)	% Go through only channels w/ significant coh
		chan = chan_sub(i);
		temp = filtfilt_rst(b, 1, emg(chan,:) ); % Band-pass EMG to frq range of interest
		temp = detrend(temp);
		bp.emg(i,:) = temp ./ (max(temp)-min(temp));	% Normalize emg data
		bp.xcor(i,:) = xcorr(bp.sdf,bp.emg(i,:),xcorr_lag,'coeff');
		pk_ind(i) = find( abs(bp.xcor(i,:)) == max( abs(bp.xcor(i,:)) ) );
		pk_val(i) = bp.xcor( i, pk_ind(i) );
		pk_lag(i) = xcorr_tlag( pk_ind(i) );
		xcor_sum(i) = sum(bp.xcor(i,:));
	end	
	
	%%%%%%%%%%%%%% 
	% Set up to plot time data
	figure
	colordef white;
	set(gcf, 'PaperPositionMode','auto','PaperOrientation','landscape');
	left = LEFT;
	width = 0.45-left;
	height = 0.90/n_xcorr;      % Used for both time & xcorr plots
	% Size to make it look good
	c = get(gcf);
	c.Position(2) = 275;
	c.Position(3) = 870;
	c.Position(4) = 680;
	set(gcf,'Position',c.Position);
	
	%%%%%%%%%%%%%% Plot normalized Spike/EMG temporal data
	left = LEFT;
	width = 0.6-left;
	height = 0.80/(n_xcorr);      % Used for both time & xcorr plots
	bottom = 0.1-height;
	% Plot EMG time domains
	for i = n_xcorr:-1:1
        bottom = bottom + height;
        subplot('position',[left bottom width height]);
		
		plot(spk.time5,bp.sdf,'Color',[0.5,0.5,0.5],'LineWidth',2);	% Plot unit data
		hold on
        plot(spk.time5,bp.emg(i,:),'k-');
        if chan_sub(i) < n_emg & i<n_xcorr
			ylabel([ 'SDF & EMG-' int2str(chan_sub(i))]);
			set(gca,'xtick',[ ]);
        elseif i == n_xcorr & chan_sub(i)~=5
			ylabel([ 'SDF & EMG-' int2str(chan_sub(i))]);
            xlabel('Time (sec)');
        elseif i == n_xcorr & chan_sub(i)==5
			ylabel('SDF & ACCEL');
            xlabel('Time (sec)');
        end
		if i == 1
			title([spk.fname ':  '  int2str(COH_FBAND(1)) '-' int2str(COH_FBAND(2)) ' Hz' ]);
		end
		ytks = get(gca,'ytick');
		ytks([1,end]) = [];
		set(gca,'ytick',ytks);
		lims = axis;
		axis([ min(spk.time5) max(spk.time5) lims(3) lims(4) ]);
		str(1) = { ['Coh =' num2str(xcor_coh(i),3)] };
		str(2) = { [ '@ ' num2str(xcor_frq(i),2) ' Hz'] };
		txt_x = lims(1) + TXT_MRG*(lims(2)-lims(1));
		txt_y = lims(4) - TXT_MRG*(lims(4)-lims(3));
		text( txt_x, txt_y, str,...
			'VerticalAlignment','top','HorizontalAlignment','left',...
			'FontSize',8);
		set(gca,'FontSize',8);
	end

	% Setup to plot xcorrs
	left = left+width+0.07;
	width = 0.95-left;
	bottom = 0.1-height;
	% Plot spike/EMG cross-correlations
	for i = n_xcorr:-1:1
        bottom = bottom + height;
        subplot('position',[left bottom width height]);
		plot(xcorr_tlag,bp.xcor(i,:),'k-');
        if i<n_xcorr
			xtks = get(gca,'xtick');
			xtks([1 2 end]) = [];
			set(gca,'xtick',xtks);
        elseif i == n_xcorr & i~=5
            xlabel('SDF Lead/Lag (msec)');
		end
		ylabel('Correlation EMG-vs-SDF');
		ytks = get(gca,'ytick');
		ytks([1,end]) = [];
		set(gca,'ytick',ytks);
		lims = axis;
		axis([ min(xcorr_tlag) max(xcorr_tlag) lims(3) lims(4) ]);
		grid on;
		str(1) = { [ 'Pk R =' num2str(pk_val(i),3) ] };
		str(2) = { [ '@ ' int2str( pk_lag(i) ) ' msec'] };
		txt_x = lims(1) + TXT_MRG*(lims(2)-lims(1));
		txt_y = lims(4) - TXT_MRG*(lims(4)-lims(3));
		text( txt_x, txt_y, str,...
			'VerticalAlignment','top','HorizontalAlignment','left',...
			'FontSize',8);
		set(gca,'FontSize',8);
	end
% 	print -dwinc
	%%%%%%%%%%%%%% END of Plotting EMG data
