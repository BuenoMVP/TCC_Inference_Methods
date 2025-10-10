function [G]=NIMCE(data,timepoint,window_num,time_lag,h,lamda2)  
[row,~]=size(data);
nsample=row/timepoint;
window_data=createTime_lag(data,nsample,window_num-1);  
% window_datatarget=createTime_lag(datatarget,nsample,window_num-1);

for ws=1:window_num
    for we=1:window_num
        if (we-ws)==time_lag
                     G_zero_order=creatG(window_data{ws},window_data{ws+1},window_data{we},window_data{ws+1},window_data{we},h);
        end
    end
end

G=edgereduce(data,G_zero_order,timepoint,window_num,time_lag,h,lamda2);

