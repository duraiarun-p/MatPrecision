clc;clear;close all;
% profile on;
% addpath('/home/s1785969/RDS/MATLAB/DAtoolbox');
% addpath('/home/s1785969/RDS/MATLAB/BlockMatchingAlgoMPEG');
% addpath('/home/s1785969/Documents/MATLAB/MotionDynamics');
% addpath('/home/arun/MATLAB/DA_Image_REG');
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
% /home/arun/Documents/MATLAB/DicomGUI/Scripting


%%
directory_name1='/home/arun/Documents/MATLAB/ImageDB/Bladder/blar9/0780334';
directory_name2='/home/arun/Documents/MATLAB/ImageDB/Bladder/blar9/D780334';
% [Cra,R2a,R1a,DSC2a,DSC1a]=BLMIR(directory_name1a,directory_name2a);
%%
m=500;
n=500;
indSize=3; % IndSize must always be odd for the symmetrical Z-directional slices
mbSize=5; % Block Size
p=5; % Search Parameter
lambda1A=0.25;% Distance Regularisation
lambda3=0.005;
costV=2^16+1;
lambda1A=lambda1A';
lambda3=lambda3'; % Orientation Regularisation
LenX=m/mbSize; % Block col index
LenY=n/mbSize; % Block row index
%%

[IP1a,IT1a,MP1a,MT1a,Z1,ZL,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
MP2ac=MP1a(TF1x);
MT2ac=MT1a(TF2x);
MP2a=MP2ac{1,1};
MT2a=MT2ac{1,1};
%%
IP1=gpuArray(double(IP1a));
IT1=gpuArray(double(IT1a));
MP1=gpuArray(MP2a);
MT1=gpuArray(MT2a);
%%
function [im1,im2,MaskV1,MaskV2,Z1,ZL,SN1,SN2]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,Imgmodal)
%Load data 1 with structure names
myhandle1=Open_DicomGui();
[~,~]=DG_load_data(myhandle1,directory_name1,Imgmodal);
[SN1,~] = DG_get_names(myhandle1);
SN1=SN1';
[im1,Z1,~]=DG_get_image_data(myhandle1);
%Load data 2 with structure names
myhandle2=Open_DicomGui();
[~,~]=DG_load_data(myhandle2,directory_name2,Imgmodal);
[SN2,~] = DG_get_names(myhandle2);
SN2=SN2';
[im2,~,~]=DG_get_image_data(myhandle2);

[~,~,ZL]=size(im1);

SN1L=length(SN1);
SN2L=length(SN2);

MaskV1=cell(SN1L,1);
for contouri=1:SN1L
    [MaskV1{contouri,1},~]=DG_generate_volume_mask(myhandle1, contouri);
end

MaskV2=cell(SN2L,1);
for contouri=1:length(SN2)
    [MaskV2{contouri,1},~]=DG_generate_volume_mask(myhandle2, contouri);
end
Close_DicomGui(myhandle1);
Close_DicomGui(myhandle2);
end