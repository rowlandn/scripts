function [pow_comp, freq] = SpikeOscil_bare( spk_t, local )
% pow = SpikeOscil_bare(spk_t)
%	A stripped-down version of SpikeOscil
%	Computes oscillatory power compensated for effects of refractory period
%
% Input:
%	spk_t - vector times of spikes (in seconds)
%	local - if defined>0, then shuffle w/in local msec [default=0]
%
% Output:
%	pow - normalized power
%	freq - matching frequency vector
%
% Created by RST, 2005-07-30
% Revised to allow local shuffling, 2005-10-10
%

global SEG_PWR
global FS       % Sampling rate
global N_SHUF	% number of ISI shuffles for control spectra

if ~exist('local','var')
	local = 0;
end
NFFT = 2^SEG_PWR;	% For spectral analysis of spike autocorr

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

isi = diff(spk.t);
isi = isi(isi>0); %make sure isi's are good
n_isi = length(isi);
rand_isi = zeros(1,n_isi);
pow_shuf=zeros( N_SHUF, len_spect);

% Calculate psd's from globaly shuffled spike trains with the same ISI
if ~local
	for i=1:N_SHUF
		rand_isi = shuffle_isi(isi);
		rand_delt = isi2delta( rand_isi );
		[pow_rand,freq_rand] = psd( rand_delt-mean(rand_delt), NFFT, FS, hann);
		pow_shuf(i,:)=pow_rand';
	end
elseif local>0
	for i=1:N_SHUF
		cum_isi=cumsum(isi);
		t_beg=0; t_end=0;
		while t_end<length(isi)
			tmp = find(cum_isi>local+50*rand);
			if (tmp)
				t_end = tmp(1);
			else
				t_end = length(isi);
			end
			rand_isi(t_beg+1:t_end) = shuffle_isi( isi(t_beg+1:t_end) );
			cum_isi=cum_isi-cum_isi(t_end);
			t_beg=t_end; 
		end
		rand_delt = isi2delta( rand_isi );
		[pow_rand,freq_rand] = psd(rand_delt-mean(rand_delt),NFFT,FS,hann);
		pow_shuf(i,:)=pow_rand';
	end
end
spk.pow_rand = mean(pow_shuf);	% Mean pow of randomized ISI's

% Now compute ISI-compensated spectrum
pow_comp = spk.pow ./ spk.pow_rand;


