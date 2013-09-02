function [f,t,cl] = sp2_fn2(d1,d2,samp_rate,flags);
% function [f,t,cl] = sp2_fn2(d1,d2,samp_rate,flags);
%
% Function with core routines for periodogram based spectral estimates.
% Implements a disjoint section analysis of two stationary signals.
% see NOTE below about calling this routine.
%
% Copyright (C) 2002, David M. Halliday.
% This file is part of NeuroSpec.
%
%    NeuroSpec is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    NeuroSpec is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with NeuroSpec; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%    NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
%
%  Inputs are two matrices containing pre processed time series or point process data.
%  Matrices have L columns with T rows.
%   L = No of disjont sections: seg_tot.
%   T = DFT segment length: seg_size.
% 
% Input arguments
%  d1         Channel 1 data matrix.
%  d2         Channel 2 data matrix.
%  samp_rate  Sampling rate (samples/sec)
%  flags      Structure with flags to control processing options.
%              Two flags supported.
%               line: Apply simple filter to suppress mains/line frequency (0:No; 1:Yes).
%               inv:  Invert channel reference for phase & cumulant        (0:No; 1:Yes).
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%
% Output parameters
%  f column 1  frequency in Hz.
%  f column 2  Log  input/d1 spectrum.
%  f column 3  Log output/d2 spectrum.
%  f column 4  Coherence.
%  f column 5  Phase.
%
%  t column 1  Lag in ms.
%  t column 2  Cumulant density.
%
%  cl.seg_size     Segment length.
%  cl.seg_tot      Number of segments.
%  cl.seg_tot_var  Effective no of segments, same as actual no for this analysis.
%  cl.samp_tot     Number of samples analysed.
%  cl.samp_rate    Sampling rate of data (samps/sec).
%  cl.dt           Time domain bin width (ms).
%  cl.df           Frequency domain bin width (Hz).
%  cl.f_c95        95% confidence limit for Log spectral estimates.
%  cl.ch_c95       95% confidence limit for coherence.
%  cl.q_c95        95% confidence limit for cumulant density.
%
% Reference:
% Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
% A framework for the analysis of mixed time series/point process data -
%  Theory and application to the study of physiological tremor,
%  single motor unit discharges and electromyograms.
% Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
%
% function [f,t,cl] = sp2_fn2(d1,d2,samp_rate,flags);
%
% NOTE: This routine is not intended to support analysis of raw data.
% It is intended as a support routine for the 2 channel spectral analysis functions:
% sp2_m.m, sp2a_m.m, sp2a2_m.m. Refer to these functions for further details.

% PBMB refers to above Progress in Biophysics article.

[seg_size,seg_tot]=size(d1); % Calculate seg_size & seg_tot from data matrices.
samp_tot=seg_tot*seg_size;   % Total no of samples, R=LT in PBMB.
fd1=fft(d1);                 % Take DFT across columns/segments ch 1, PBMB (4.1)/(4.2).
fd2=fft(d2);                 % Take DFT across columns/segments ch 2, PBMB (4.1)/(4.2). 

t_fac=2*pi*seg_size;               % Normalization for periodogram spectral estimates.
                                   % NB 1/L included in mean() function below.
f11=mean(abs(fd1.*fd1)/t_fac,2);   % Spectrum 1, PBMB (5.2), Mag squared for auto spectra.
f22=mean(abs(fd2.*fd2)/t_fac,2);   % Spectrum 2, PBMB (5.2), Mag squared for auto spectra.
f21=mean(fd2.*conj(fd1)/t_fac,2);  % Cross spectrum (complex valued),PBMB (5.2).
deltaf=samp_rate/seg_size;         % Resolution - spacing of Fourier frequencies in Hz.

% Set line frequency in Hz.
%line_freq=50;
line_freq=60;  % Uncomment this line to set the line frequency to 60 Hz.

% Suppression of mains/line frequency - smooth out using adjacent values.
if flags.line      
  line_ind=round(line_freq/deltaf)+1;        % NB Index 1 is DC.
  f11(line_ind)=0.5*(f11(line_ind-2)+f11(line_ind+2));    % Spectrum ch 1.
  f11(line_ind-1)=0.5*(f11(line_ind-2)+f11(line_ind-3));
  f11(line_ind+1)=0.5*(f11(line_ind+2)+f11(line_ind+3));
  f22(line_ind)=0.5*(f22(line_ind-2)+f22(line_ind+2));    % Spectrum ch 2.
  f22(line_ind-1)=0.5*(f22(line_ind-2)+f22(line_ind-3));
  f22(line_ind+1)=0.5*(f22(line_ind+2)+f22(line_ind+3));
  f21(line_ind)=0.5*(f21(line_ind-2)+f21(line_ind+2));    % Cross spectrum.
  f21(line_ind-1)=0.5*(f21(line_ind-2)+f21(line_ind-3));
  f21(line_ind+1)=0.5*(f21(line_ind+2)+f21(line_ind+3));
% Smooth elements in upper hermetian section of cross spectral estimate.
% This data used in ifft() to generate cumulant. Data is complex conjugate.
  f21(seg_size-line_ind+2)=conj(f21(line_ind));
  f21(seg_size-line_ind+3)=conj(f21(line_ind-1));
  f21(seg_size-line_ind+1)=conj(f21(line_ind+1));
end

% Channel reference invert: Channel 2 is now reference(input) for phase & cumulant
%  Spectra and coherence unaffected by this.
if flags.inv
  f21=conj(f21);
end 

% Construct output spectral matrix f.
seg_size_2=(2:seg_size/2+1)';    % Indexing for output, DC component not output.
f(:,1)=(seg_size_2-1)*deltaf;    % Column 1 - frequencies in Hz.
f(:,2)=log10(f11(seg_size_2));   % Column 2 - Log spectrum ch 1.
f(:,3)=log10(f22(seg_size_2));   % Column 3 - Log spectrum ch 2.
                                 % Column 4 - Coherence, PBMB (5.5).
f(:,4)=abs(f21(seg_size_2)).*abs(f21(seg_size_2))./(f11(seg_size_2).*f22(seg_size_2));
f(:,5)=angle(f21(seg_size_2));   % Column 5 - Phase, PBMB (5.7).

% Estimate cumulant density using inverse DFT of cross spectrum.
deltat=1000.0/samp_rate;    % dt in msec.

cov=ifft(f21);              % Inverse DFT.

% Construct output time domain matrix t.
% Column 1 - time in msec. Range (-T/2)*dt to (T/2-1)*dt.
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;

% Column 2 - cumulant, shifted by T/2 so that time zero is in centre.
% 2pi/T factor is 2*pi, since ifft routine includes 1/T term.
t([seg_size/2+1:seg_size,1:seg_size/2],2)=real(cov(1:seg_size))*2*pi; % PBMB (5.9).

% Estimate variance of cumulant density estimate.
var_fac=4*pi*pi/(seg_size*samp_tot);                           % Factor (2pi/T)(2pi/R).
q_var=var_fac*2*sum(f11(1:seg_size/2+1).*f22(1:seg_size/2+1)); % PBMB (6.10).

% Construct cl structure, confidence limits for parameter estimates.
cl.seg_size=seg_size;    % T.
cl.seg_tot=seg_tot;      % L.
cl.seg_tot_var=seg_tot;  % Effective no of segments: for multivariate & pooled spectra.
cl.samp_tot=samp_tot;    % R.
cl.samp_rate=samp_rate;  % Sampling rate.
cl.dt=deltat;            % Delta t.
cl.df=deltaf;            % Delta f.
cl.f_c95=0.8512*sqrt(1/seg_tot);    % 95% Confidence limit for spectral estimates, PBMB (6.2).
% N.B. Confidence interval for log plot of spectra is TWICE this value.
cl.ch_c95=1-0.05^(1/(seg_tot-1));   % 95% Confidence limit for coherence, PBMB (6.6).
cl.ch_c999=1-0.001^(1/(seg_tot-1));   % 99.98% Confidence limit for coherence, PBMB (6.6).
cl.q_c95=1.96*sqrt(q_var);          % 95% Confidence limits for cumulant, PBMB (6.11).
