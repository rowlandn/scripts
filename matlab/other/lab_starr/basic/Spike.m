load gaus_flt
spk_len = length(ad0);
flt_len = (length(gaus_50ms_20kHz)-1)/2;
spk = (diff(ad0 < -0.07)>0);
spk(spk_len)=0;
smspk = conv(gaus_50ms_20kHz,spk);
smspk(1:flt_len)=[];
smspk(spk_len+1:spk_len+flt_len)=[];