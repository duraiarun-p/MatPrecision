clc;clear;close all;
%%
rsiz=100;    
Rm=randi(20,rsiz);
Rm1=randi(200,rsiz);
%%
N=hist3([Rm(:),Rm1(:)],'CdataMode','auto');
xlabel('Rm')
ylabel('Rm1')
colorbar
view(2)