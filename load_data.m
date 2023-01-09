% load target data of the 4 recording sites 
%    ref{s}.CSD:     CSD   [200*5 timepoints x 12 channels]
%    ref{s}.MUA:     MUA   [200*5 timepoints x 16 channels]
%    ref{s}.times:   times [1 x 200 timepoints]
%    ref{s}.freqs:   5 tones (BF tone and 4 non-BF tones) 
%    ref{s}.delayIn: input delay (msec)
%    ref{s}.MUApeak: peak value of MUA 

function ref=load_data()
    path='.\data\';
    ref=cell(4,1);
    tmp1=load([path, 'site1.mat']);
    tmp2=load([path, 'site2.mat']);
    tmp3=load([path, 'site3.mat']);
    tmp4=load([path, 'site4.mat']);
    ref{1}=tmp1.ref;
    ref{2}=tmp2.ref;
    ref{3}=tmp3.ref;
    ref{4}=tmp4.ref;
end