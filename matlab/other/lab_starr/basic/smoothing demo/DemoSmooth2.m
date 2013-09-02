function DemoSmooth2(n,h)
% Re-draws graph when number of passes slider is moved
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

Passes=round(n);
DemoSmoothRedraw