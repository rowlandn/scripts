function SmoothSlider1(n,h)
% Re-draws graph when smooth width slider is moved
% Tom O'Haver, July 2006. % Slider function by Matthew Jones.
global X
global SmoothY
global Y
global SmoothWidth
global Passes
SmoothWidth=round(n);
SmoothSliderRedraw