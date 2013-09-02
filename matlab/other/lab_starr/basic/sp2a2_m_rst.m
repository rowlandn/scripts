function [f,t,cl] = sp2a2_m(dat1,dat2,samp_rate,seg_pwr,opt_str);
% function [f,t,cl] = sp2a2_m(dat1,dat2,samp_rate,seg_pwr,opt_str)
%
% Function to calculate spectra, coherence, phase & cumulant for 2 time series
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
%  dat1       Channel 1  (input) time series array (multiple columns, rst mod).
%  dat2       Channel 2 (output) time series array (multiple columns, rst mod).
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
%  r  Rectification - this options requires an argument, valid arguments are:[0, 1, 2].
%                      r0:  Rectify Input  channel  (dat1)
%                      r1:  Rectify Output channel  (dat2)
%                      r2:  Rectify Both   channels (dat1 & dat2)
%  t  linear De-trend - this options requires an argument, valid arguments are:[0, 1, 2].
%                      t0:  De-trend Input  channel  (dat1)
%                      t1:  De-trend Output channel  (dat2)
%                      t2:  De-trend Both   channels (dat1 & dat2)
% Options examples:
%  to set all options on, set opt_str='i m r2 t2'
%  to rectify and de-trend both channels, set opt_str='r2 t2' 
%  to rectify channel 2 and de-trend both channels, set opt_str='r1 t2'
%  to rectify channel 1 and de-trend both channels, set opt_str='r0 t2'
%
% Output parameters
%  f column 1  frequency in Hz.
%  f column 2  Log input/dat1  spectrum.
%  f column 3  Log output/dat2 spectrum.
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
%  cl.opt_str      Copy of options string.
%  cl.what         Text label, used in plotting routines.
%
%
%	RST modification:  allow dat1 & dat2 to be arrays (multiple columns) of
%	continuously-sampled data, disjoint between columns
%
% Reference:
% Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
% A framework for the analysis of mixed time series/point process data -
%  Theory and application to the study of physiological tremor,
%  single motor unit discharges and electromyograms.
% Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
%
% function [f,t,cl] = sp2a2_m(dat1,dat2,samp_rate,seg_pwr,opt_str)

% Check numbers of arguments.
if (nargin<4)
  error(' Not enough input arguments');
end  

if (nargout<3)
  error(' Not enough output arguments');
end  

% Check for single column data
[nrow1,ncol1]=size(dat1);
[nrow2,ncol2]=size(dat2);
if (ncol1~=ncol2 | nrow1~=nrow2)
  error('Data arrays are not equal size')
end 

pts_tot= nrow1*ncol1;           % Determine size of data vector.
if ((nrow2*ncol2)~=pts_tot)      % Check that input vectors are equal length. 
  error (' Unequal length data arrays');
end
seg_size=2^seg_pwr;             % DFT segment length (T).
seg_col=fix(nrow1/seg_size);	% Number of segments per column
seg_tot=ncol1*seg_col;  		% Total Number of complete segments (L).
samp_tot=seg_col*seg_size;      % Number of samples to analyse per column.

% Issue warning if number of segments outside reasonable range,
%  may indicate a problem with analysis parameters.
if (seg_tot<6)
  warning(['You have a small number of segments: ',num2str(seg_tot)])
end 
if (seg_tot>1000)
  warning(['You have a large number of segments: ',num2str(seg_tot)])
end 

% Arrange data into L columns each with T rows.
rd1=reshape(dat1(1:samp_tot,:),seg_size,seg_tot);
rd2=reshape(dat2(1:samp_tot,:),seg_size,seg_tot);
md1=mean(rd1);      % Determine mean of each column/segment.
md2=mean(rd2);

% Process options.
flags.line=0;       % Set defaults - options off.
flags.inv=0;
trend_chan_1=0; 
trend_chan_2=0;
rect_chan_1=0;
rect_chan_2=0;
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
    case 'r'             % Rectification option.        
      i=str2num(optarg);
    if (i<0 | i>2)
      error(['error in option argument -- r',optarg]);
    end  
    if (i~=1)
      rect_chan_1=1;     % Rectify ch 1.
    end  
    if (i>=1)
      rect_chan_2=1;     % Rectify ch 2.
    end  
    case 't'             % Linear de-trend option.
      i=str2num(optarg);
    if (i<0 | i>2)
      error(['error in option argument -- t',optarg]);
    end  
    if (i~=1)
      trend_chan_1=1;    % De-trend ch 1.
    end  
    if (i>=1)
      trend_chan_2=1;    % De-trend ch 2.
    end  
    case 'm'             % Mains/line frequency suppression option.
      flags.line=1;
    otherwise
      error (['Illegal option -- ',opt]);  % Illegal option.
  end
end

if (trend_chan_1 | trend_chan_2)       % Index for fitting data with polynomial.
  trend_x=(1:seg_size)';
end

for ind=1:seg_tot                        % Loop across columns/segments.
  rd1(:,ind)=rd1(:,ind)-md1(ind);          % Subtract mean from ch 1.
  rd2(:,ind)=rd2(:,ind)-md2(ind);          % Subtract mean from ch 2.
  if rect_chan_1
    rd1(:,ind)=abs(rd1(:,ind));            % Rectification of ch 1 (Full wave).
  end
  if rect_chan_2
    rd2(:,ind)=abs(rd2(:,ind));            % Rectification of ch 2 (Full wave).  
  end
  if trend_chan_1                               % Linear trend removal.
    p=polyfit(trend_x,rd1(:,ind),1);            % Fit 1st order polynomial.
    rd1(:,ind)=rd1(:,ind)-p(1)*trend_x(:)-p(2); % Subtract from ch 1.
  end  
  if trend_chan_2                               % Linear trend removal.
    p=polyfit(trend_x,rd2(:,ind),1);            % Fit 1st order polynomial.
    rd2(:,ind)=rd2(:,ind)-p(1)*trend_x(:)-p(2); % Subtract from ch 2.
  end  
end

% Call sp2_fn2() Periodogram based spectral estimation routine.
[f,t,cl]=sp2_fn2(rd1,rd2,samp_rate,flags);

% Set additional elements in cl structure.
cl.N1=0;                              % N1, No of events in ch 1. (zero for TS data)
cl.N2=0;                              % N2, No of events in ch 2.          "
cl.P1=0;                              % P1, mean intensity ch1.            "
cl.P2=0;                              % P2, mean intensity ch2.            "
cl.opt_str=opt_str;                   % Copy of options string.
cl.what='';                           % Field for plot label.

% Display No of segments & resolution 
disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),' sec,  Resolution: ',num2str(cl.df),' Hz.']);
