function DemoSmooth1(n,h)
% Re-draws graph when smooth width slider is moved
% Tom O'Haver, July 2006. Slider function by Matthew Jones.

global t
global SmoothSignal
global signal
global Noise
global NoiseArray
global puresignal
global wid
global SmoothWidth
global Passes

SmoothWidth=round(n);
DemoSmoothRedraw