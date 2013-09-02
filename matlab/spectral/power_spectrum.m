function mean_psd = power_spectrum(data,nfft,Fs,sp1,sp2,sp3)

% power_spectrum This function performs a power spectrum operation
% on given data.  Data can be either single or multiple
% traces in column vector format.  The number of fast-fourier 
% transforms (nfft) should be a power of 2 (e.g., 512,1024, etc) and
% should be equal to the number of data points that would 
% encompass at least 5 cycles of oscillations in the data trace(s).
% The traces are resampled to 400Hz using the original sampling rate (Fs) 
% of the recordings in Hz.
% 
% mean_psd = power_spectrum(data,nfft,Fs);
%
% Example 1: mean_psd = power_spectrum(EEG_L,1024,10000);

for i = 1:size(data,2)
    
    [all_psd(:,i) F] = psd(data(:,i),nfft,Fs,'mean');%*(1/Fs);
    
end

mean_psd.mean_psd = mean(all_psd,2);
mean_psd.freq = F;

%figure
if nargin == 6
    subplot(sp1,sp2,sp3)
else
end
%plot(mean_psd.freq,mean_psd.mean_psd)
semilogy(mean_psd.freq,mean_psd.mean_psd)
