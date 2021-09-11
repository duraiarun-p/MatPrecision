% clc;clear;close all;
%%
a=1:100;
b=0:100:9900;
A=zeros(100);
%%
for i=1:100
    A(i,:)=b(i)+a;
end