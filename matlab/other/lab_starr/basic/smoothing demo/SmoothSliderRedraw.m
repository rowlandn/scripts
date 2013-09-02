% Redraws graph for SmoothSlider when the sliders are changed
  axes(h);
  PlotRange=[SmoothWidth.*1:length(X)-SmoothWidth.*1];
  temp=Y;
  for k=1:Passes,
      % You can use any smooth function here in place of fastsmooth,
      temp=fastsmooth(temp, SmoothWidth);
  end
  SmoothY=temp;
  plot(X(PlotRange),SmoothY(PlotRange))
  figure(1);
  title(['Smooth Width = ' num2str(SmoothWidth) '   Number of passes = ' num2str(Passes) '    Signal maximum = ' num2str(max(SmoothY)) ])
  xlabel('The smoothed signal is in the vector SmoothY')
  h2=gca;
  axis([X(1) X(length(X)) min(Y) max(Y)]);
  grid on