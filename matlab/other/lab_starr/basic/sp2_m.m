function [f,t,cl] = sp2_m(sp1,sp2,samp_rate,seg_pwr,opt_str);
% function [f,t,cl] = sp2_m(sp1,sp2,samp_rate,seg_pwr,opt_str)
%
% Function to calculate spectra, coherence, phase & cumulant for 2 spike trains
% using periodogram based estimation procedure for stationary signals.
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
% Input arguments
%  sp1        Channel 1   (input) spike times - specify as integer no of sampling intervals,
%                                               ordered in ascending order.
%  sp2        Channel 2  (output) spike times - specify as integer no of sampling intervals,
%                                               ordered in ascending order.
%  samp_rate  Sampling rate (samples/sec).
%  seg_pwr    Segment length - specified as power of 2.
%  opt_str    Options string.
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%
% Input Options
%  i  Invert channel reference for phase and cumulant density.
%  m  Mains/line frequency suppression.
% Options examples:
%  to set both options on, set opt_str='i m'
%
% Output parameters
%  f column 1  frequency in Hz.
%  f column 2  Log  input/sp1 spectrum.
%  f column 3  Log output/sp2 spectrum.
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
%  cl.N1           N1, No of events/spikes in ch1. 
%  cl.N2           N2, No of events/spikes in ch2.
%  cl.P1           P1, mean intensity ch1.
%  cl.P2           P2, mean intensity ch2.
%  cl.opt_str      Copy of options string.
%  cl.what         Text label, used in plotting routines.
%
% Reference:
% Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
% A framework for the analysis of mixed time series/point process data -
%  Theory and application to the study of physiological tremor,
%  single motor unit discharges and electromyograms.
% Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
%
% function [f,t,cl] = sp2_m(sp1,sp2,samp_rate,seg_pwr,opt_str)

% Check numbers of arguments
if (nargin<4)
  error(' Not enough input arguments');
end  

if (nargout<3)
  error(' Not enough output arguments');
end

% Check for single column data
[nrow,ncol]=size(sp1);
if (ncol~=1)
  error(' Input NOT single column: sp1')
end 
[nrow,ncol]=size(sp2);
if (ncol~=1)
  error(' Input NOT single column: sp2')
end 

% Check spike data for any negative values.
if (min(sp1)<1 | min(sp2)<1)
  error(' Negative or zero spike times in input data');
end

% Check spike data for any non integer values.
if ~isequal(round(sp1),sp1)
  error(' Non integer spike times in sp1');
end
if ~isequal(round(sp2),sp2)
  error(' Non integer spike times in sp2');
end

% Size of data array - taken as longest spike train.
pts_tot=max([sp1(length(sp1)) sp2(length(sp2))]);

% Set up 0/1 point process representations.
dat1(1:pts_tot,1)=0;
dat1(sp1)=1;
dat2(1:pts_tot,1)=0;
dat2(sp2)=1;

seg_size=2^seg_pwr;             % DFT segment length (T)
seg_tot=fix(pts_tot/seg_size);  % Number of complete segments (L)
samp_tot=seg_tot*seg_size;      % Number of samples to analyse: R=LT.

% Issue warning if number of segments outside reasonable range,
%  may indicate a problem with analysis parameters.
if (seg_tot<6)
  warning(['You have a small number of segments: ',num2str(seg_tot)])
end 
if (seg_tot>1000)
  warning(['You have a large number of segments: ',num2str(seg_tot)])
end 

% Arrange data into L columns each with T rows.
rd1=reshape(dat1(1:samp_tot),seg_size,seg_tot);
rd2=reshape(dat2(1:samp_tot),seg_size,seg_tot);
md1=mean(rd1);      % Determine mean of each column/segment.
md2=mean(rd2);

% Process options.
flags.line=0;       % Set defaults - options off.
flags.inv=0;
if (nargin<5)  % No options supplied.
  opt_str='';
end
options=deblank(opt_str);
while (any(options))              % Parse individual options from string.
  [opt,options]=strtok(options);
  optarg=opt(2:length(opt));      % Determine option argument.
  switch (opt(1))
    case 'i'             % Channel reference invert option.
      flags.inv=1;
    case 'm'             % Mains/line frequency suppression option.
      flags.line=1;
    otherwise
      error (['Illegal option -- ',opt]);  % Illegal option.
  end
end

for ind=1:seg_tot                         % Loop across columns/segments.
  rd1(:,ind)=rd1(:,ind)-md1(ind);          % Subtract mean from ch 1.
  rd2(:,ind)=rd2(:,ind)-md2(ind);          % Subtract mean from ch 2.
end

% Call sp2_fn2() Periodogram based spectral estimation routine.
[f,t,cl]=sp2_fn2(rd1,rd2,samp_rate,flags);

% Set additional elements in cl structure.
cl.N1=length(sp1);                    % N1, No of events in ch 1. 
cl.N2=length(sp2);                    % N2, No of events in ch 2.
cl.P1=cl.N1/(pts_tot*2*pi);           % P1, mean intensity ch1.
cl.P2=cl.N2/(pts_tot*2*pi);           % P2, mean intensity ch2.
cl.opt_str=opt_str;                   % Copy of options string.
cl.what='';                           % Field for plot label.

% Display No of segments & resolution 
disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),' sec,  Resolution: ',num2str(cl.df),' Hz.']);
