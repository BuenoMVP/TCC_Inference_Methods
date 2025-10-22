
%data=data;
clear; clc;
%rime clock
tStart = tic;
comparedata={'net3_expression_data';};

% Load TSV data
data = dlmread('/home/marco/projects/TCC_Inference_Methods/Database/input data/net3_expression_data.tsv', '\t');
% datatarget=colom_634_target_normal_dataexp;
window_num=3;
time_lag=2;
h=0.45;
timepoint=21;
[G]=NIMCE(data,timepoint,window_num,time_lag,h,0.04);
dlmwrite(['/home/marco/projects/TCC_Inference_Methods/NIMCE-master/database','NIMCE_',comparedata{1},'.txt'], G, 'delimiter', '\t','precision',100,'newline', 'pc')
% running time
toc(tStart)
