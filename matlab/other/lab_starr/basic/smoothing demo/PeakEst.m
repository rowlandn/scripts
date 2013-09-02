function [PeakX, PeakY, Width]=PeakEst(x,y,pos,wid)
% Least-squares estimate of peak height, position, and width of a Gaussian
% peak in the signal x,y located at approximately x=pos and with approximate 
% half-width "wid".  Fits parabola to log of points near the peak.
% [PeakX, PeakY, Width]=PeakEst(x,y,pos,wid)
% Tom O'Haver, July 2006
   xrange=[(pos-wid):(pos+wid)];
   coef=polyfit(x(xrange),log(abs(y(xrange))),2);  % Fit parabola to log of sub-group
   c1=coef(3);c2=coef(2);c3=coef(1);
   PeakX=-c2/(2*c3);   % Compute peak position and height of Gaussian
   PeakY=exp(c1-c3*(c2/(2*c3))^2);
   Width=2.35703/(sqrt(2)*sqrt(-1*c3));