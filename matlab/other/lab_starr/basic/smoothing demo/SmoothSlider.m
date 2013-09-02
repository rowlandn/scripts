% Interactive smoothing, with slider to control Smooth Width
% Place the signal to be smoothed in the global variables X,Y, and 
% define MaxSmoothwidth as the largest value of the smooth width
% slider. Smoothed signal is left in global variable SmoothY.
% The actual smoothing is performed by the function SmoothSliderRedraw,
% which is called when the sliders are moved. You can change 
% the smoothing function used by editing SmoothSliderRedraw.  
% Tom O'Haver, toh@umd.edu, July 2006.  Slider function by Matthew Jones.
figure(1);
close

global X
global SmoothY
global Y
global MaxSmoothwidth
global SmoothWidth
global Passes
SmoothWidth=1;
Passes=1;
  h=figure(1);
  plot(X,Y)
  figure(1);title(['Original Signal.   Signal maximum = ' num2str(max(Y)) ])  
  grid on
  h2=gca;axis([X(1) X(length(X)) min(Y) max(Y)]);
  rtslid(h,@SmoothSlider1,h2,1,'Scale',[1 MaxSmoothwidth],'Def',1,'Back',[0.9 0.9 0.9],'Label','Smooth');
  rtslid(h,@SmoothSlider2,h2,1,'Scale',[1 10],'Def',1,'Back',[0.9 0.9 0.9],'Label','Passes','Position',[0.95 0.1 0.03 0.8]);