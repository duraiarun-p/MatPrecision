%% Motion Estimation - Exhaustive Search based Block Matching - Choose Second Method
clc;clear;close all;
% profile on;
% addpath('/home/s1785969/RDS/MATLAB/DAtoolbox');
% addpath('/home/s1785969/RDS/MATLAB/BlockMatchingAlgoMPEG');
% addpath('/home/s1785969/Documents/MATLAB/MotionDynamics');
addpath('/home/arun/Documents/MATLAB/MAWS_Precision');
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
% /home/arun/Documents/MATLAB/DicomGUI/Scripting


%%

directory_name1='/home/arun/Documents/MATLAB/ImageDB/Bladder/blar1/0428705';
directory_name2='/home/arun/Documents/MATLAB/ImageDB/Bladder/blar1/D428705';

%% Main Module to call BLM function for mismatch data & needed corrections

[IP1a,IT1a,MP1a,MT1a,Z1a,~,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
MP2ac=MP1a(TF1x);
MT2ac=MT1a(TF2x);
MP2a=MP2ac{1,1};
MT2a=MT2ac{1,1};
%%
ind=33;
sigma=2;
I1=IP1a(:,:,ind);
I2=IT1a(:,:,ind);
I1s=imgaussfilt(I1,sigma);
I2s=imgaussfilt(I2,sigma);
Ie=imabsdiff(I1,I2);
Ies=imabsdiff(I1s,I2s);
[Gx1,Gy1]=gradient(double(I1s));
[Gx2,Gy2]=gradient(double(I2s));
Gex=imabsdiff(Gx1,Gx2);
Gey=imabsdiff(Gy1,Gy2);
Ge=Gex+Gey;
% C=normxcorr2(I1,I2);
% CGx=normxcorr2(Gx1,Gx2);
% CGy=normxcorr2(Gy1,Gy2);
% CG=CGx+CGy;
figure(10);subplot(121),imshow(I1,[]);subplot(122),imshow(I1s,[]);
figure(20);subplot(121),imshow(I2,[]);subplot(122),imshow(I2s,[]);
figure(30),subplot(121),imshow(Ie,[]);subplot(122),imshow(Ies,[]);
figure(40),subplot(121),imshow(Gex,[]);subplot(122),imshow(Gey,[]);
figure(50),subplot(121),imshow(Ies,[]);title('Smooth');subplot(122),imshow(Ge,[]);title('Summed Gradient');
% figure;subplot(131),imshow(Gx1,[]);subplot(132),imshow(Gx2,[]);subplot(133),imshow(CGx,[]);
% figure;subplot(131),imshow(Gy1,[]);subplot(132),imshow(Gy2,[]);subplot(133),imshow(CGy,[]);
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
end
