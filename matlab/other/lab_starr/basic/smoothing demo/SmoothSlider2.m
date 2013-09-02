function SmoothSlider2(n,h)
% Re-draws graph when smooth passes slider is moved
% Tom O'Haver, July 2006. % Slider function by Matthew Jones.
global X
global SmoothY
global Y
global SmoothWidth
global Passes
Passes=round(n);
SmoothSliderRedraw