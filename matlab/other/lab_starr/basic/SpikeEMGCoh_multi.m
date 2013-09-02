function SpikeEMGCoh_multi()
% SpikeEMGCoh_multi()
% This program analyses Spike times versus EMG data
% Created by RST, 2002-04-05
% Modified by RST, 2003-10-15
% Modified by RST, 2005-07-12 to handle multiple spikes in each NEX file
%
%	Run within a directory and this program will process all spikes in all 
%	nex files in the directory

clear all

global FS       % Sampling rate for processed EMG & spike SDF
global FSS		% Subsample SDF & EMG to this rate
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
global COH_FBAND
global XCORR_LAG

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

spkname_pattern = 'Snip\w*[abcd]';


cd(uigetdir);

FileLst = dir('*.nex');
if(length(FileLst) < 1)
	error(['Found no NEX files in current directory.' ]);
end

% Open output file and print header line
	fname = 'SpikeEMGCoh.txt';
	outfile = fopen(fname,'w');
	if(outfile == -1)
       error(['Unable to open...' fname ]);
	end
	fprintf(outfile,'fname\trecord_dur\t');
	for i = 1:N_ACorr
		% For max reported signif acorr pks, header freq, SNR & oscill index
		fprintf(outfile,'AC%d_frq\tAC%d_SNR\tAC%d_OI\t',i,i,i);
	end
	fprintf(outfile,'AC_more?\tXCoh_SIG\t');	% Header for missed significant ACorr peaks & XCoh 95% CI
	for i = 1:N_XCorr
		fprintf(outfile,'XCoh%d_chan\tXCoh%d_frq\tXCoh%d_coh\tXCoh%d_phs\t',i,i,i,i);
	end
	fprintf(outfile,'XCoh_more?\t');		% Header for missed significant XCoh peaks
	for i = 1:N_XCorr
		fprintf(outfile,'XCor%d_chan\tXCor%d_cor\tXCor%d_lag\tXCor%d_sum\t',i,i,i,i);
	end
	fprintf(outfile,'\n');		% EOL

	for m = length(FileLst):-1:1
		if findstr(FileLst(m).name,'_m.nex')	% Unable to process merged files
			FileLst(m) = [];
		end
	end

m = 0;
while( m < length(FileLst))
	m = m+1;
	disp( ['Processing... ' FileLst(m).name] );
	
    [nvar, varname, types] = nex_info(FileLst(m).name);
	
	[spk{m}, emg{m}, emg_cl{m}] = SpikeEMG_multi(FileLst(m).name,...
		spkname_pattern);
	
%	print(gcf);
	

	fprintf(outfile,'%s\t%.3f\t',FileLst(m).name,spk{m}.dur);
	% Report significant acor peaks
	for i = 1:N_ACorr
		if i <= length(spk{m}.sig_pks)
				
			fprintf(outfile,'%.3f\t%.3f\t%.3f\t',...
				spk{m}.frq(spk{m}.sig_pks(i)),...
				spk{m}.snr(spk{m}.sig_pk_ind(i)),...
				spk{m}.oscil(spk{m}.sig_pk_ind(i)) );
		else
			fprintf(outfile,'\t\t\t');
		end
	end
	if length(spk{m}.sig_pks) > N_ACorr
		fprintf(outfile,'yes\t');
	else
		fprintf(outfile,'no\t');
	end	
	
	if ~isnan(emg{m})
		[d,d,n_emg] = size(emg{m});
		% Find significant xcor peaks
		frq_ind = find( emg{m}(:,FRQ_COL,1) < SRCH_HI );
		xcor_thresh = emg_cl{m}.coh_sig;
		fprintf(outfile,'%.3f\t',xcor_thresh);
		coh = emg{m}(frq_ind,COH_COL,:);
		phs = emg{m}(frq_ind,PHS_COL,:);
		frq = emg{m}(frq_ind,FRQ_COL,:);

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
			fprintf(outfile,'yes\t');
		else
			fprintf(outfile,'no\t');
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%
		% Do XCorr for sig coh peak w/in COH_FBAND
		xinds = find(xcor_frq > COH_FBAND(1) & xcor_frq < COH_FBAND(2));
		n_xcors = length(xinds);
		if n_xcors
			[pk_val, pk_lag, xcor_sum] = SpikeEMGXCorr(FileLst(m).name,chan_sub(xinds),xcor_coh(xinds),xcor_frq(xinds));

			for i = 1:N_XCorr
				if i <= n_xcors
					fprintf(outfile,'%d\t%.3f\t%.3f\t%.3f\t',...
						chan_sub(xinds(i)),pk_val(i),pk_lag(i),xcor_sum(i));
				else
					fprintf(outfile,'\t\t\t\t');
				end
			end
			fprintf(outfile,'\n');
		else
			for i = 1:N_XCorr
				fprintf(outfile,'\t\t\t\t');
			end
			fprintf(outfile,'\n');
		end		
	else
		fprintf(outfile,'\t');
		for i = 1:N_XCorr
			fprintf(outfile,'\t\t\t\t');
		end
		fprintf(outfile,'\t');
		for i = 1:N_XCorr
			fprintf(outfile,'\t\t\t\t');
		end
		fprintf(outfile,'\n');
	end
		
end
fclose(outfile);

save SpikeEMGCoh
