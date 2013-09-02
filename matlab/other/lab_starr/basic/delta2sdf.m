function sdf = delta2sdf(delta,gaus_sigma)
% sdf = delta2sdf(delta,gaus_sigma)
%
% Convert delta array to mean spike density function (sdf)
%  Argument:	delta = trial X msec array,  1 denoting spike events
%				gaus_sigma = width of gaussian (in msec; optional, default = 20)
%
%  Returns:		sdf = mean spike density function across trials
%
%	04-01-27	Revised to take delta in either column or row vector
%	Revised 2005-08-24 to make gaussian on the fly - any sigma can be used
%			also speeded execution by first summing delta and convolving
%			one summed-delta stream.

[ntrial spk_len] = size(delta);

% Make the filter
if ~exist('gaus_sigma','var')
	gaus_sigma = 20;
end
flt_len = 5*gaus_sigma;	% usually 5*stdev is long enough to make is smooth
x = -1*flt_len:flt_len;
flt = normpdf(x,0,gaus_sigma);
% Make sure edges are smooth & scaling correct (weights should sum to 1000)
flt = flt - min(flt);
flt = 1000.*flt ./(sum(flt)*ntrial);

flt_len = (length(flt)-1)/2;

siz = size(delta);

if min(siz) == 1
	ntrial = 1;
	spk_len = max(siz);
	if siz(2) == 1
		delta = delta';
	end
else
	ntrial = siz(1);
	spk_len = siz(2);
end

% sum the delta function across trials
sdelta = sum(delta,1);

% indices into right and left ends of original delta vector
revlist1 = flt_len+1:-1:2;
revlist2 = spk_len-1:-1:spk_len-flt_len;

% build a mirror-padded sdelta
s = [ sdelta(revlist1) sdelta sdelta(revlist2) ];

% convolve w/ filter
sdf = conv(flt, s );

% remove pads
sdf(1: 2*flt_len)=[];
sdf(spk_len+1: spk_len + 2*flt_len)=[];

if siz(2) == 1
	sdf = sdf';
end
