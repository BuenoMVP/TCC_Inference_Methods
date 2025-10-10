pkg load statistics;

%data=data;
clear; clc;
%rime clock
tStart = tic;
comparedata={'insilico_size100_1_timeseries'};

%ȡ
load (comparedata{1});
% load (comparedata{2});
data=insilico_size100_1_timeseries; 
% datatarget=colom_634_target_normal_dataexp;
%㷨
window_num=3;
time_lag=2;
h=0.45;
timepoint=21;
[G]=NIMCE(data,timepoint,window_num,time_lag,h,0.04);
%Ϊһԡ\tΪtxtļע·
dlmwrite(['NIMCE_',comparedata{1},'.txt'], G, 'delimiter', '\t','precision',100,'newline', 'pc')
% running time
toc(tStart)