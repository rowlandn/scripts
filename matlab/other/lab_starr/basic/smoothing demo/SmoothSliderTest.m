% Test of SmoothSlider
% Generates a synthetic signal assigned to Y, then calls SmoothSlider.
% Tom O'Haver, toh@umd.edu, July 2006

global X
global SmoothY
global Y
global MaxSmoothwidth
global SmoothWidth
global Passes

X=[1:10000]; 
Noise=1;
Y=3.*gaussian(X,3000,1500)+gaussian(X,7000,1000)+Noise.*randn(size(X)); % Generate synthetic test signal
MaxSmoothwidth=500;  % Defines the largest value of the smooth width slider.
SmoothSlider