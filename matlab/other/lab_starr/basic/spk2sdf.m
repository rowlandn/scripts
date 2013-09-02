function [SDF] = spk2sdf(spk)
% [SDF] = spk_to_sdf(spk)
%
% Convert spk array to mean spike density function (SDF)
%  Argument:	spk = trial X msec array,  1 denoting spike events
%
%  Returns:		SDF = mean spike density function across trials

load gaus_flt
siz = size(spk);
ntrial = siz(1);
spk_len = siz(2);
flt_len = (length(gaus_10ms_1kHz)-1)/2;

% Pad spike vectors at ends to avoid boundary smoothing to zero
s = zeros(1, spk_len + 2*flt_len);
% indices into right and left ends of original spk vector
revlist1 = flt_len+1:-1:2;
revlist2 = spk_len-1:-1:spk_len-flt_len;
SDF = zeros(1,spk_len);
for i = 1:ntrial,
	s = [ spk(i,revlist1) spk(i,:) spk(i,revlist2) ];
	smspk = conv(gaus_20ms_1kHz, s );
	smspk(1: 2*flt_len)=[];
	smspk(spk_len+1: spk_len + 2*flt_len)=[];
	SDF = SDF + smspk;
end

SDF = SDF/ntrial;
