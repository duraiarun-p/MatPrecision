clc;clear;close all;
addpath('/home/arun/Documents/MATLAB/ImageOutputs/PerfComp');
%%
load('Perf_1.mat','T','subfolder');
%% Dont change it works for now

Sco=cell(1,5);
for fj=1:length(subfolder)
Sc=zeros(1,4);
for fi=2:length(subfolder{fj,1})
    Pn{fi,1}=subfolder{fj,1}(fi).name;
    s=T{fj,fi};
    Sc=[Sc;s];
end
Sc(1,:)=[];
Pn(1,:)=[];
Pn=[Pn';Pn'];
Pn=reshape(Pn,[numel(Pn) 1]);
Sc=num2cell(Sc);
Sc=[Pn,Sc];
Sco=[Sco;Sc];
clear Sc Pn s
end
Sco(1,:)=[];
Tab=cell2table(Sco);
writetable(Tab,'Score.xlsx');
%% parameter value can only be changed in original file for easier analysis
% Sco=cell(1,5);
% fj=1;
% Sc=zeros(1,4);
% for fi=2:length(subfolder{fj,1})
%     Pn{fi,1}=subfolder{fj,1}(fi).name;
%     s=T{fj,fi};
%     Sc=[Sc;s];
% end
% Sc(1,:)=[];
% Pn(1,:)=[];
% Pn=[Pn';Pn'];
% Pn=reshape(Pn,[numel(Pn) 1]);
% Sc=num2cell(Sc);
% Sc=[Pn,Sc];
% Sco=[Sco;Sc];
% % clear Sc Pn s
% 
% Sco(1,:)=[];
% Tab=cell2table(Sco);
% writetable(Tab,'Score.xlsx');