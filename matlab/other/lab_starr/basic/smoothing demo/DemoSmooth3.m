function DemoSmooth3(n,h)
% Re-draws graph when the Resample slider is moved
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

NoiseArray=Noise.*randn(size(t)); % Generate another random noise samplpe
signal=puresignal+NoiseArray; % Add noise to peak
DemoSmoothRedraw