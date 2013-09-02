% this script plots the frequency spectrum of files containing digitized
% ecog data.  The script looks for paired variables named "ecogrest" and
% "ecogactive" (recordings from same contact montage, during rest and active self-cued movement)
% and plots the time domain and frequency domain for each, side by side

Fs=1000; %FS = samples per unit time, in this case digitization at 1000 hz rate, may pull from glr file
Gain=5700; %amplifier gain on ecog channels,will modify program to make the gain an input, or pull gain from 
% glr file
z=length(ecogrest);
lastnum=(z-1)/1000;
ecogresttime=0:0.001:lastnum;

subplot(2,2,1);
plot(ecogresttime,ecogrest/Gain); %plots the resting ecog in time domain, NOTE not yet clear what voltage 
%scale is after correcting for gain, see email to anthony
xlabel('time (seconds)');
ylabel('voltage');
title('ecog at rest in time domain')

Y=fft(ecogrest);
Y(1)=0; %for some reason the first number in the fft vector is the sum of all others, 
% and to avoid confusion in the plotting, we set that first number to zero
n=length(Y); 
P=(abs(Y).^2)/n; %calculates the power of the fft, which is square of its magnitude
P=P./max(P); %normalizes the power to the maximum power of 1 to facilate plotting
P1=P(1:floor(n./2));%pulls out only the first half of the spectrum up to nyquist frequency
n1=length(P1);
f=(0:n1-1)*(Fs/n); %sets the frequency or x scale to be correct range to correspond to fft
subplot(2,2,3);
plot(f,P1); %plots the power versus the frequency
axis([0 100 0 1]); %expands the first 100 Hz of the spectrum
xlabel('frequency (Hz)');
ylabel('normalized power');
title('ecog frequencies at rest, 0-100 hz')


za=length(ecogactive);
lastnuma=(za-1)/1000;
ecogactivetime=0:0.001:lastnuma;

subplot(2,2,2);
plot(ecogactivetime,ecogactive/Gain); %plots the active movement ecog in time domain
xlabel('time (seconds)');
ylabel('voltage');
title('ecog during movement in time domain')

Ya=fft(ecogactive);
Ya(1)=0; %for some reason the first number in the fft vector is the sum of all others, 
% and to avoid confusion in the plotting, we set that first number to zero
na=length(Ya); 
Pa=(abs(Ya).^2)/na; %calculates the power of the fft, which is square of its magnitude
Pa=Pa./max(Pa); %normalizes the power to the maximum power of 1 to facilate plotting
P1a=Pa(1:floor(na./2));%pulls out only the first half of the spectrum up to nyquist frequency
n1a=length(P1a);
fa=(0:n1a-1)*(Fs/na); %sets the frequency or x scale to be correct range to correspond to fft
subplot(2,2,4);
plot(fa,P1a); %plots the power versus the frequency
axis([0 100 0 1]); %expands the first 100 Hz of the spectrum
xlabel('frequency (Hz)');
ylabel('normalized power');
title('ecog frequency during movement, 0-100 hz')




