clc;clear;close all;
%%
a=1:2020;a=a';
batch_size=20;
num_of_batch=round(length(a)/batch_size);

for index = 1:num_of_batch
%     start_index(index)=index*batch_size;
%     stop_index(index)=min((index+1)*batch_size,length(a));
stop_index(index)=min(index*batch_size,length(a));
start_index(index)=stop_index(index)-batch_size+1;
end