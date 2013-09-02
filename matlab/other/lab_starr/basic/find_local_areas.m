function [area, pk_ind] = find_local_areas(Y,thresh,peaks)
% [area, pk_ind] = find_local_areas(Y,thresh,peaks)
% This program finds the area above thresh for each peak in vector X
% and returns indices of their peak positions
%
%	Usage:  [area, pk_ind] = find_local_areas(Y,thresh,peaks)
%		where 	Y = vector of arbitrary length defining curve
%				thresh = compute area under curve & above threshold
%				area = areas of each peak
%				pk_ind = index of peak value for each area
% Created by RST, 2003-10-15
% 			Revised 2004-01-27
xpand = 10;		% expand input vector by this factor to get more accurate

YY = interp1(Y,(1/xpand):(1/xpand):length(Y))';

YY = YY - thresh;

peaks = peaks .* xpand;	% Indices into expanded YY vector
npeaks = length(peaks);
% 1st find start indices for each peak
for i = 1:npeaks
	if i==1
		sep = find( YY(1:peaks(i)) < 0 );
		if isempty(sep)
			start(i) = 1;
		else
			start(i) = max(sep)+1;
		end
	else
		sep = find( YY(peaks(i-1):peaks(i)) < 0 ) + peaks(i-1) - 1;
		if isempty(sep)
			start(i) = find(  YY(peaks(i-1):peaks(i)) == min(YY(peaks(i-1):peaks(i))) ) + peaks(i-1) - 1;
		else
			start(i) = max(sep)+1;
		end
	end
end

% 1st - find start indices for each peak
for i = 1:npeaks
	if YY( peaks(i) ) <= 0 
		start(i) = peaks(i);
		continue;
	elseif i==1
		sep = find( YY(1:peaks(i)) < 0 );
		if isempty(sep)
			start(i) = 1;
		else
			start(i) = max(sep)+1;
		end
	else
		sep = find( YY(peaks(i-1):peaks(i)) < 0 ) + peaks(i-1);
		if isempty(sep)
			start(i) = find( YY(peaks(i-1):peaks(i)) == min(YY(peaks(i-1):peaks(i))) ) + peaks(i-1) - 1;
		else
			if max(sep)<peaks(i)
				start(i) = max(sep)+1;
			else
				start(i) = peaks(i);
			end
		end
	end
end

% 2nd - find stop indices for each peak
for i = 1:npeaks
	if YY( peaks(i) ) <= 0 
		stop(i) = peaks(i);
		continue;
	elseif i==npeaks
		sep = find( YY(peaks(i):length(YY)) < 0 ) + peaks(i) - 1;
		if isempty(sep)
			stop(i) = length(YY);
		else
			stop(i) = min(sep)-1;
		end
	else
		sep = find( YY(peaks(i):peaks(i+1)) < 0 ) + peaks(i) - 1;
		if isempty(sep)
			stop(i) = find( YY(peaks(i):peaks(i+1)) == min(YY(peaks(i):peaks(i+1))) ) + peaks(i) - 1;
		else
			if min(sep)>peaks(i)
				stop(i) = min(sep)-1;
			else
				stop(i) = peaks(i);
			end
		end
	end
end

% 3rd - Compute area under each peak
area = zeros(1,npeaks);
for i = 1:npeaks
	if start(i)==stop(i)
		continue;
	else
		area(i) = sum( YY( start(i):stop(i) ) );
	end
end
area = area ./ xpand;

return