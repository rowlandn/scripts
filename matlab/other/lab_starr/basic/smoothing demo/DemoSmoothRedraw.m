% Redraws graph for DemoSmoothSlider when the sliders are changedglobal t
global SmoothSignal
global signal
global Noise
global NoiseArray
global puresignal
global wid
global SmoothWidth
global Passes
  axes(h);
  PlotRange=[SmoothWidth.*1:length(t)-SmoothWidth.*1];
  
 % Smooth the simulated signal with noise
  temp=signal;
  for k=1:Passes,
      % You can use any smooth function here in place of fastsmooth
      temp=fastsmooth(temp, SmoothWidth);
  end
  SmoothSignal=temp;
  
 % Smooth the noise
  temp=NoiseArray;
  for k=1:Passes,
      % You can use any smooth function here in place of fastsmooth
      temp=fastsmooth(temp, SmoothWidth);
  end
  SmoothNoise=temp;
  
 % Smooth the pure noiseless signal
  temp=puresignal;
  for k=1:Passes,
      % You can use any smooth function here in place of fastsmooth
      temp=fastsmooth(temp, SmoothWidth);
  end
  PureSmoothSignal=temp;
  % Compute the SNR
  Signal=range(PureSmoothSignal(PlotRange));
  SNR=Signal./std(SmoothNoise(PlotRange));  
  % Plot the smoothed signal
  h=figure(1);
  plot(t(PlotRange),SmoothSignal(PlotRange))
  % Compute the peak parameters of the smoothed peak
  [PeakX, PeakY, Width]=PeakEst(t,SmoothSignal,length(t)./2,wid/2);
  figure(1);title(['Smooth Width = ' num2str(SmoothWidth) '     Smooth Ratio = ' num2str(SmoothWidth./wid) '     Passes = ' num2str(Passes) ])
  xlabel(['Signal maximum = ' num2str(max(SmoothSignal))  '     Peak Width = ' num2str(Width) '     SNR = ' num2str(SNR) ]);
  h2=gca;
  axis([0 length(t) -.10 1]);
  grid on