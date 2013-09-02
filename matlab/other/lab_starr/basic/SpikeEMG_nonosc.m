function [result] = SpikeEMG_nonosc(nexname)
% [result] = SpikeEMG_nonosc(fname)
% This program analyses Spike times versus EMG data
% Created by RST, 2003-11-02
%
%	
clear all
global FS       % Sampling rate for processed EMG & spike SDF
corr_leadlag = 5000;	% Look for correlations up to this lead & lag (msec relative to 0=synchronous)
BOOT = 1000;
FS = 1000;

low_passes = [ 0.2, 0.5, 1, 2, 5 ];	% Compute correlations on data low-passed @ these freq's

color_lst = [ 'r', 'g', 'b', 'm', 'c', 'k'];

deci_rates = round( FS ./ low_passes );

[nexname, RTPATH] = uigetfile('*.nex', 'Select spike time file (NEX)');
if (nexname ~= 0)
    [nvar, varname, types] = nex_info(nexname);
    if nvar > 1
        error([num2str(nvar) ' spikes in ' nexname '!!!  I don''t know how to process > 1 spike yet']);
    end
    [spk.n, spk.t] = nex_ts(nexname,varname);
else
    error(['I can''t find the NEX file:  ' nexname ]);
end

fname = strrep(nexname,'.nex','');

tmp = dir([fname '.mat']);
[sz x] = size(tmp);
if(sz ==0)
    error(['I can''t find the EMG file:  ' tmp ]);
else
	EMG_FLAG = true;
	load( [fname '.mat'] );
	[n_emg, data_len] = size(emg_chan);
end
emg_chan(5,:) = abs(emg_chan(5,:));

% Convert spike times to delta function & sdf
spk.delt = spk_t2delta(spk.t,data_len);
spk.sdf = spk2sdf(spk.delt);

% 
for i = 1:length(deci_rates)
	d_spk = decimate( spk.sdf, deci_rates(i) );
	result{i}.spk = d_spk;	
	for j = 1:n_emg
		d_emg = decimate( emg_chan(j,:), deci_rates(i) );
		result{i}.emg(j,:) = d_emg;	

		% Calc correlations at different lead/lags
		leadlag_n = round(corr_leadlag ./ deci_rates(i));
		home_ind = (leadlag_n+1):length(d_spk)-leadlag_n;
		for k = -leadlag_n:leadlag_n
			
			[R,P] = corrcoef(d_spk(home_ind),d_emg(home_ind+k));
			result{i}.corr(j,k+leadlag_n+1) = R(1,2);
			result{i}.pear(j,k+leadlag_n+1) = P(1,2);
			result{i}.phs(j,k+leadlag_n+1) = k*deci_rates(i)./FS;
		
			if exist('bootstrp')
				% Bootstrap to better estimate significance of correlation
				rhos = bootstrp(BOOT,'corrcoef',d_spk(home_ind),d_emg(home_ind+k));
				rhos = sort(rhos(:,2));
				if R(1,2) < 0
					tail = find( rhos >= 0 );
					if isempty(tail)
						result{i}.bootp(j,k+leadlag_n+1) = 0;
					else
						result{i}.bootp(j,k+leadlag_n+1) = 1 - min(tail) ./ BOOT;
					end
				else
					tail = find( rhos <= 0 );
					if isempty(tail)
						result{i}.bootp(j,k+leadlag_n+1) = 0;
					else
						result{i}.bootp(j,k+leadlag_n+1) = max(tail) ./ BOOT;
					end
				end
			else
				result{i}.bootp(j,k+leadlag_n+1) = NaN;
			end
		end
	end
end

figure
for j = 1:n_emg
	subplot(n_emg,1,j);
	hold on;
	for i = 1:length(deci_rates)
		plot(result{i}.phs(j,:),result{i}.corr(j,:),color_lst(i))
		set(gca,'FontSize',8);
	end
	str1 = sprintf('Spk-vs-EMG%d (R)',j);
	ylabel(str1);
	xlabel('EMG lead/lag (sec)');
	if j == 1
		title(nexname,'FontSize',12);
	end
end
for i = 1:length(deci_rates)
	str{i} = sprintf('%.2f Hz',low_passes(i));
end
legh = legend(str);
set(legh,'FontSize',8);
