clc;clear;close all;
% profile on;
addpath('/home/arun/Documents/MATLAB/MAWS_Precision');
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');

mypath='/home/arun/Documents/MATLAB/ImageDB/PrintoutDB/';
outpath='/home/arun/Documents/MATLAB/ImageOutputs/DicomOutputs';
%%
patfolder=dir(mypath);
patfolder(1:2,:) = []; 
bigfolderlen=length(patfolder);
start=1;
stop=start;
% stop=2;
% stop=bigfolderlen;
plotrow=round(sqrt(stop));
if (plotrow^2<stop)
    plotcol=plotrow+1;
else
    plotcol=plotrow;
end
tic;
for folderi=start:stop
mypath1=[patfolder(folderi).folder,'/',patfolder(folderi).name];
flag=0;
if flag==0
load(mypath1,'ContourName','PlanCTMaskcell','PlanCTPath','PlanCTInfo','PlanCTData','PlanCTLoca','PreCBCTInfo','PreCBCTPath','PreCBCTData','PreCBCTLoca');
else
load(mypath1,'ContourName','PlanCTMaskcell','PlanCTPath','PlanCTInfo','PlanCTData','PlanCTLoca','PosCBCTInfo','PosCBCTPath','PosCBCTData','PosCBCTLoca');
end
% Space for batch processing
end
timeL=toc;
%% Selecting bladder mask
maskind1=find(strcmpi(ContourName,'Bladder'));
PlanMask=PlanCTMaskcell{maskind1,1};
%% Rigid Registration stage
tic;
[geomtform,Rfixed,Rmoving,optimizer,metric]=dblmrigidim(PreCBCTData,PlanCTData,PreCBCTInfo{1,1},PlanCTInfo{1,1});
timeR=toc;
PlanCTR = imwarp(PlanCTData,Rmoving,geomtform,'bicubic','OutputView',Rfixed);%Warping scan 
PlanMask1 = imwarp(PlanMask,Rmoving,geomtform,'bicubic','OutputView',Rfixed);%Warping mask
%%
figure(1),
subplot(131),
imshow(PreCBCTData(:,:,round(length(PreCBCTLoca)/2)),[])
title('Pre CBCT scan slice @ mid-plane')
subplot(132),
imshow(PlanCTData(:,:,round(length(PlanCTLoca)/2)),[])
title('Plan CT scan slice @ mid plane');
subplot(133),
imshow(PlanCTR(:,:,round(length(PreCBCTLoca)/2)),[])
title('After rigid registration')
%%
figure(2),
subplot(131),
histogram(PreCBCTData);
title('Histogram-Pre CBCT scan')
subplot(132),
histogram(PlanCTData)
title('Histogram-Plan CT scan');
subplot(133),
histogram(PlanCTR)
title('Histogram-Rigidly registered Plan CT scan')
% subplot(235),
% imshow(PlanMask(:,:,round(length(PlanCTLoca)*0.15)),[])
% title('Plan CT mask @ mid plane');
% subplot(236),
% imshow(PlanMask1(:,:,round(length(PreCBCTLoca)*0.15)),[])
% title('After rigid registration')