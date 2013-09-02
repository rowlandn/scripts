function psp_ch1(f,cl,freq,ch_max,label)
% function to plot coherence in current figure/subplot window
%  psp_ch1(f,cl,freq,ch_max,label)
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
%		 Contact:  contact@neurospec.org
%
%  psp_ch1(f,cl,freq,ch_max,label)
%
% f,cl     Output from spectral analysis routine.
% freq     Frequency limit for plotting (Hz).
% ch_max   Maximum value of y axis (Optional).
% label    Optional title instead of cl.what.
global MAX_N

frq_col = 1;	% coherence in 4th column
coh_col = 4;	% coherence in 4th column
phs_col = 5;	% phse in 5th column
 

frq_step = (f(2,frq_col)-f(1,frq_col))/2;

freq_pts=round(freq/cl.df);
phase_pts = zeros(freq_pts*2,1);
phase_frq = zeros(freq_pts*2,1);
phase_pts(:) = NaN;
phase_frq(:) = NaN;
coh_pks = find_local_nmax( f(1:freq_pts,coh_col), MAX_N);
sig_pts = find(f(1:freq_pts,coh_col) >= cl.coh_sig);
sig_pks = intersect( coh_pks, sig_pts);
phase_label = [];
for i = 1:length(sig_pts)
	phase_pts(sig_pts(i)*2-1) = f(sig_pts(i),phs_col);
	phase_pts(sig_pts(i)*2) = f(sig_pts(i),phs_col);
	phase_frq(sig_pts(i)*2-1) = f(sig_pts(i),frq_col)-frq_step;
	phase_frq(sig_pts(i)*2) = f(sig_pts(i),frq_col)+frq_step;
	if ismember( sig_pts(i), sig_pks)
		phase_label = [phase_label, f(sig_pts(i),phs_col)];
	end
end
ind = find(isfinite(phase_pts));
if length(ind)
	phaseplt_mx = max(phase_pts(ind))+20;
	phaseplt_mn = min(phase_pts(ind))-20;
else
	phaseplt_mx = 1;
	phaseplt_mn = -1;
end

%Check freq range
[x,y]=size(f);
if (freq_pts>x)
	error('Requested frequency range too large.');
end

f_max=freq_pts*cl.df;
ch_cl=[cl.df,cl.coh_sig;f_max,cl.coh_sig];

h1 = plot(f(1:freq_pts,frq_col),f(1:freq_pts,coh_col),'k',...
	ch_cl(:,1),ch_cl(:,2),'k:','LineWidth',2);
ax1 = gca;
set(ax1,'Box','off');
if (nargin<4)
  axis([0,freq,-Inf,Inf]);
else
  axis([0,freq,0,ch_max]);
end
xtks = get(ax1,'XTick');
xtks(:,1) = [];
xtks(:,end) = [];
set(ax1,'XTick',xtks);
ytks = get(ax1,'YTick');
ytks(:,1) = [];
set(ax1,'YTick',ytks,'FontSize',8);
xlabel('Freq (Hz)');
ylabel('R^2');

ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom',...
	'YAxisLocation','right','Color','none','YColor','r');
	
% Plot phase on second axis
h2 = line(phase_frq,phase_pts,...
	'Color','r','LineWidth',0.3,'Parent',ax2);
axis([0,freq,phaseplt_mn,phaseplt_mx]);
set(ax2,'XTick',[]);
ytks = sort( unique(round(phase_label)) )';
dif_tk = diff( ytks);
if min( dif_tk ) < 5
	ytks( find(dif_tk == min(dif_tk))+1) = [];
end	
set(ax2,'YTick',ytks,'FontSize',8);
ylabel('\phi (ms, EMG-SDF)');
if (nargin>4)
  title(label);
else
 title(['Coherence: ',cl.what]);
end  
