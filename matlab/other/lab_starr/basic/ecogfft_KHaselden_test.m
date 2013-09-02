

Fs=1000;            % samples per unit time, in this case digitization at 1000 hz rate

% fft_ecog_signal performs FFT on raw ecog signal and outputs normalized power
% in the frequency domain up to the Nyquist frequency

Y=fft(raw_ecog_signal);
% will use syntax Y=fft(ecogsnippet1, N) where N is the number of points that are
% used for the time domain (and frequency domain), use N=2000)
% note in this section we will take the starting times for each rest
% snippet or each active snippet, run the fft for each.
Y(1)=0; %for some reason the first number in the fft vector is the sum of all others,
% and to avoid confusion in the plotting, we set that first number to zero
n=length(Y);
P=(abs(Y).^2)/n; %calculates the power of the fft, which is square of its magnitude
P_out = P(1:floor(n/2));
% at this point, do the quantitative analyses - many of the functions below
% this point (normalizing to peak power, smoothing, taking logs) are really
% just to facilitate plotting and data visualization

%lets try the following for each individual channel, in each individual movement type (jaw, hand etc)
%The five "raw" P vectors for each epoch should be averaged together (prior to
%any normalization).  then "normalize" the movement related averaged P vector by the averaged rest
%P vector. Plot normalized P vector analagous to fig 2 of liu and aziz 2008) Then segment the normalized P vector into frequency bins
%(bandwidth to be determined, but could start with 10 hz bands 0 to 150)
%and determine mean normalized power in that frequency bin.  write these
%data to a table as in the example.  
%rest vector 

norm_P=P./max(P); %normalizes the power to the maximum power of 1 to facilate plotting
fft_out = norm_P(1:floor(n/2));%pulls out only the first half of the spectrum up to nyquist frequency
% n1 = length(fft_out);
% freq_out =(0:n1-1)*(Fs/n); %sets the frequency or x scale to be correct range to correspond to fft
% % run the above sequence for each of the 5-10 snippets at rest, and each in
% % movement, and add the fft_out vectors for each to a final
freq_unsmoothed = linspace(0,Fs/2,length(fft_out));
plot(freq_unsmoothed, fft_out);
xlim([0 50]);
ylim([0 1]);
title('FFT RCalvert snapshot 1 ecog #1 - #6');