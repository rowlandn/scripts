% Demonstration of smoothing on a single peak, with sliders to control
% smooth width and the number of passes of the smoothing function through the
% signal. Generates a synthetic signal with random noise, smooths it, and
% measures the effect on signal-to-noise ratio (SNR), peak height, and peak 
% width. The Resample slider takes another random noise sample. To change 
% the soothing function, replace "fastsmooth" in SmoothSampleRedraw.
% Tom O'Haver, toh@umd.edu, July 2006. Slider function by Matthew Jones.
figure(1);
close
global t
global SmoothSignal
global signal
global Noise
global NoiseArray
global puresignal
global wid
global SmoothWidth
global Passes

t=[1:2000];
amp=1 ;    % Amplitude of the peak
pos=1000;  % Position of the peak
wid=200;   % Width of the peak
Noise=.3;
% Generate peak signal with noise
puresignal=amp.*gaussian(t,pos,wid); 
NoiseArray=Noise.*randn(size(t)); % Generate random noise
signal=puresignal+NoiseArray; % Add noise to peak
% Compute initial signal-to-noise ratio
Signal=range(puresignal);
SNR=Signal./std(NoiseArray);
SmoothWidth=1;
Passes=1;
h=figure(1);
% Plot the simulated signal
plot(t,signal)
[PeakX, PeakY, Width]=PeakEst(t,signal,length(t)./2,wid/2);
figure(1);title(['Smooth Width = ' num2str(SmoothWidth) '     Smooth Ratio = ' num2str(SmoothWidth./wid) '     Passes = ' num2str(Passes) ])
xlabel(['Signal maximum = ' num2str(max(signal))  '     Peak Width = ' num2str(Width) '     SNR = ' num2str(SNR) ]);
grid on
h2=gca;axis([0 length(t) -.1 1]);
% Draw the sliders
MaxSmoothwidth=200; % Maximukl range of SmoothWidth slider (change if desired)
rtslid(h,@DemoSmooth1,h2,1,'Scale',[1 MaxSmoothwidth],'Def',1,'Back',[0.9 0.9 0.9],'Label','Smooth');
rtslid(h,@DemoSmooth2,h2,1,'Scale',[1 10],'Def',1,'Back',[0.9 0.9 0.9],'Label','Passes','Position',[0.95 0.5 0.03 0.35]);
rtslid(h,@DemoSmooth3,h2,1,'Scale',[0 1],'Def',0,'Back',[0.9 0.9 0.9],'Label','Resample','Position',[0.95 0.05 0.03 0.35]);
 