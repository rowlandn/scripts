function [avg, st_err] = quantCOH(A, FREQ_QCOH, freq)

% matrix A is an aggregated coherence data with data from each patient
% stored in its columns

% initialize variables
% coh has 5 rows for 5 frequency bands, and a column for each patient
coh = zeros(length(FREQ_QCOH),size(A,2));
avg = zeros(size(FREQ_QCOH,1),1);
st_err = zeros(size(FREQ_QCOH,1),1);

for i=1:size(FREQ_QCOH,1)
    idx = find(freq>FREQ_QCOH(i,1) & freq<FREQ_QCOH(i,2));
    coh(i,:) = mean(A(idx,:),1); 
    avg(i) = mean(coh(i,:));
    st_err(i) = std(coh(i,:))/sqrt(size(A,2));
end