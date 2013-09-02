%  Make Test data for SpikeEMGCoh

clear all

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

SEG_PWR = 9;	% Segment length - specified as power of 2
FS = 1000;
SRCH_LO = 2.0;	% low freq for search for significant oscillations in autocorr
SRCH_HI = 30;	% top freq for search for significant oscillations
SIG_SNR = 7;	% significant oscill must have power SNR > SIG_SNR * std(power)  :7 => false pos in 5% of cases
SIG_OSC = 10;	% or OSC_IND > SIG_OSC  :10 => false pos in 5% of cases
SIG_COH = 0.0002;	% Significance threshold for Spike-EMG coherence :0.0002 => false positive in 5% of cases

MAX_N = 7;		% find max's in acorr over MAX_N symmetric points (must be odd)
AUTOCORR_LAG = 500;

% Indices into columns of Spike/EMG data (returned from NeuroSpec)
FRQ_COL = 1;	% frequency in 1st column
COH_COL = 4;	% coherence in 4th column
PHS_COL = 5;	% phase in 5th column

N_ACorr = 2;	% max # of significant autocorr pks reported in txt output
N_XCorr = 5;	% max # of significant spike/EMG coh pks reported in txt output

PLOTTING = false;
TEST_REPS = 10000;

%%% Parameters for generating synthetic data
FRQ1 = 4;
MAG1 = 0.0;	% Percent modulation of spikes & EMG
FRQ2 = 10;
MAG2 = 0.0;
MEAN_RATE = 50;	% average spike rate of delta fn
NOISE = 1;
NEMG = 5;
EMG_PHASE_SHIFT = 0.02;

FS = 1000;
DUR = 60;	% 60 second of data


% Load filter information
load('abfconv_filters');

% Open output file and print header line
	fname = 'SpikeEMG_test4.txt';
	outfile = fopen(fname,'w');
	if(outfile == -1)
       error(['Unable to open...' fname ]);
	end
	fprintf(outfile,'fname\trecord_dur\t');
	for i = 1:N_ACorr
		% For max reported signif acorr pks, header freq, SNR & oscill index
		fprintf(outfile,'AC%d_frq\tAC%d_SNR\tAC%d_OI\t',i,i,i);
	end
	fprintf(outfile,'AC_more?\tXC_99.9CI\t');		% Header for missed significant XCorr peaks & XCorr 99.9% CI
	for i = 1:N_XCorr
		fprintf(outfile,'XC%d_chan\tXC%d_frq\tXC%d_coh\tXC%d_phs\t',i,i,i,i);
	end
	fprintf(outfile,'XC_more?\n');		% Header for missed significant XCorr peaks


m = 0;
while( m < TEST_REPS)
	m = m+1;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Make test data
	disp( ['Making test #... ' int2str(m) ] );

	time = (1/FS:1/FS:DUR);
	npts = length(time);
	
	oscil = 1 + MAG1.*sin(2*pi*FRQ1*time) + MAG2.*sin(2*pi*FRQ2*time);
	
	rand('state',sum(100*clock));
	ran_spk =  rand(1,npts);		% uniform random 0-to-1
	ran_spk = ran_spk .* oscil ;  % Delta fn w/ mean rate = 500 Hz
	
	thresh = 1 -  MEAN_RATE ./ 1000 ;
	
	spk.delt = zeros(1,npts);
	spk_ind = find( ran_spk > thresh );
	spk.delt(spk_ind) = 1;
	spk.t = find(spk.delt) .* 1/FS;		% Convert to spike times
	
	randn('state',sum(100*clock));
	
	emg_noise = NOISE .* randn(NEMG,npts) ./ 4;		% EMG white noise scaled 
	emg_oscil = zeros(NEMG,npts);		% 
	
	% NOTE:  emg_chan 1 lags spk oscil, emg_chan 5 leads spk oscil
	for i = 1:NEMG
		t = time - (NEMG-1)*EMG_PHASE_SHIFT/2 + (i-1)*EMG_PHASE_SHIFT;
		emg_oscil(i,:) = 1 + MAG1.*sin(2*pi*FRQ1*t) + MAG2.*sin(2*pi*FRQ2*t);
		emg_chan(i,:) = filtfilt( LowPass_100Hz_1kHz.tf.num, 1, abs(emg_noise(i,:) .* emg_oscil(i,:)) );
	end
	save SpikeEMG_test	
	clear spk
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Process test data
	disp( ['Processing test #... ' int2str(m) ] );
	[spk, emg, emg_cl] = SpikeEMG5_test('SpikeEMG_test.mat');
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Save test results
	disp( ['Saving test #... ' int2str(m) ] );
	fprintf(outfile,'%s\t%.3f\t',[ 'Test#' int2str(m)], spk.dur);
	% Report significant acor peaks
	for i = 1:N_ACorr
		if i <= length(spk.sig_pks)
				
			fprintf(outfile,'%.3f\t%.3f\t%.3f\t',...
				spk.frq(spk.sig_pks(i)),...
				spk.snr(spk.sig_pk_ind(i)),...
				spk.oscil(spk.sig_pk_ind(i)) );
		else
			fprintf(outfile,'\t\t\t');
		end
	end
	if length(spk.sig_pks) > N_ACorr
		fprintf(outfile,'yes\t');
	else
		fprintf(outfile,'no\t');
	end	
	
	if ~isnan(emg)
		[d,d,n_emg] = size(emg);
		% Find significant xcor peaks
		frq_ind = find( emg(:,FRQ_COL,1) < SRCH_HI );
		xcor_thresh = emg_cl.coh_sig;
		fprintf(outfile,'%.3f\t',xcor_thresh);
		coh = emg(frq_ind,COH_COL,:);
		phs = emg(frq_ind,PHS_COL,:);
		frq = emg(frq_ind,FRQ_COL,:);

		xcor_ind = find( coh > xcor_thresh);	% Find coh > threshold across all emg chans
		
		% Find local max's in coh for all emg chans (expressed as indices into
		% whole array) 
		d_len = length(coh(:,:,1));
		pk_inds = find_local_nmax(coh(:,:,1),MAX_N);
		for i = 2:n_emg
			pk_inds = [ pk_inds, (i-1)*d_len + find_local_nmax(coh(:,:,i),MAX_N) ];
		end
		
		% Find overlap
		xcor_ind = intersect(xcor_ind,pk_inds);
		
		[frq_sub,chan_sub] = ind2sub(size(coh),xcor_ind);
		xcor_coh = coh(xcor_ind);
		xcor_phs = phs(xcor_ind);
		xcor_frq = frq(xcor_ind);
	
		% Print significant xcor peaks
		for i = 1:N_XCorr
			if i <= length(xcor_coh)
				fprintf(outfile,'%d\t%.3f\t%.3f\t%.3f\t',...
					chan_sub(i),xcor_frq(i),xcor_coh(i),xcor_phs(i));
			else
				fprintf(outfile,'\t\t\t\t');
			end
		end
		if length(xcor_coh) > N_XCorr
			fprintf(outfile,'yes\n');
		else
			fprintf(outfile,'no\n');
		end
	else
		fprintf(outfile,'\t');
		for i = 1:N_XCorr
			fprintf(outfile,'\t\t\t\t');
		end
		fprintf(outfile,'\n');
	end
	clear spk;		
end
fclose(outfile);
