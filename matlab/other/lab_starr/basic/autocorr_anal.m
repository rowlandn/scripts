function [snr, oscil, pow, frq, pk_ind, acorr_ss, lag_ss] = autocorr_anal(acorr, lag, fs)
% function [snr, oscil, pow, frq] = autocorr_anal(acorr, lag, fs)
%
% Find peaks in autocorr spectrum & determin SNR & OscilIndex for each
% Inputs -
%	acorr	- autocorrelogram, symmetric around zero, msec bins
%	lag		- matching lag times for acorr
%	fs		- sampling frequency of autocorr
%
% Outputs -
% 	snr		- signal-to-noise for each spectral peak
% 	oscil	- oscillation index for each spectral peak
% 	pow		- power of autocorr spectrum
% 	frq		- matching frequency vector
%	pkind	- indices into pow & frq for spectral peaks
%	acorr_ss - sub-sampled autocorrelation
%	lag_ss	- matching sub-sampled lag
%
%	RST 05/07/18

% Load filter information
load('abfconv_filters');

% Defines from calling function - Determines significance thresholds
global SEG_PWR	
global SRCH_LO
global SRCH_HI
global SIG_SNR
global SIG_OSC	% Oscillation index thresh for significance
global MAX_N	% Max search range in autocorr local max search

if ~exist('SEG_PWR','var')
	SEG_PWR = 9;
end
if ~exist('SRCH_LO','var')
	SRCH_LO = 2.0;
end
if ~exist('SRCH_HI','var')
	SRCH_HI = 30;
end
if ~exist('MAX_N','var')
	MAX_N = 9;
end

SRCH_LO = 2.0;
SRCH_HI = 30;
MAX_N = 9;

SS_INTERV = 1000/fs;
NFFT = 2^SEG_PWR;	% For spectral analysis of spike autocorr
AUTOCORR_LAG = max(lag);

% find & remove central valley
avg_autocorr = mean(acorr(1:(AUTOCORR_LAG-10)));
for i = 1:AUTOCORR_LAG
	if acorr(i+AUTOCORR_LAG+1) >= avg_autocorr
		break
	end
end
for j = -i:i
	acorr(j+AUTOCORR_LAG+1) = avg_autocorr;
end
acorr = detrend(acorr);
temp = filtfilt(LowPass_100Hz_1kHz.tf.num, 1, acorr ); % Chop out > 100 Hz 
acorr_ss = decim(temp,SS_INTERV);	% sub-sample autocorr down to 200 Hz
lag_ss = decim(lag,SS_INTERV);	% sub-sample lags down to 200 Hz

% Do spectral computation @ 200 Hz 
[Px, frq, units, Sxx] = periodogram(acorr_ss,[],NFFT, fs);
pow = Sxx( 1:length(Px) );
acorr(1:AUTOCORR_LAG/SS_INTERV) = [];
lag(1:AUTOCORR_LAG/SS_INTERV) = [];

% Find peaks in power spectrum 
srch_ind = find( frq > SRCH_LO & frq < SRCH_HI );
pow_pks = find_local_nmax( pow, MAX_N);
[pk_ind, ind_pow_pks, ind_srch_ind] = intersect( pow_pks, srch_ind );

% Compute SNR for each peak
snr = (pow(pk_ind) - mean( pow(srch_ind))) ./ std(pow) ;

% Compute Oscil index for each peak
pow_area = find_local_areas( pow, mean( pow(srch_ind)), pow_pks );
oscil = 100 .* pow_area(ind_pow_pks)' ./ sum(pow);		% Oscillation index
