clc;clear;close all;
% profile on;

addpath('/home/arun/Documents/MATLAB/MAWS_Precision');
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
%%
mypath='/home/arun/Documents/MATLAB/ImageDB/PrintoutDB/';
patfolder=dir(mypath);
patfolder(1:2,:) = []; 
bigfolderlen=length(patfolder);

start=10;
stop=start;
% stop=2;
% stop=bigfolderlen;

plotrow=round(sqrt(stop));
if (plotrow^2<stop)
    plotcol=plotrow+1;
else
    plotcol=plotrow;
end
figure(1),
% figtitle('3D CBCT vs CT volume location match test'),
for folderi=start:stop
mypath1=[patfolder(folderi).folder,'/',patfolder(folderi).name];
flag=0;
if flag==0
load(mypath1,'ContourName','PlanCTMaskcell','PlanCTPath','PlanCTInfo','PlanCTData','PlanCTLoca','PreCBCTInfo','PreCBCTPath','PreCBCTData','PreCBCTLoca');
else
load(mypath1,'ContourName','PlanCTMaskcell','PlanCTPath','PlanCTInfo','PlanCTData','PlanCTLoca','PosCBCTInfo','PosCBCTPath','PosCBCTData','PosCBCTLoca');
end
% Space for batch processing
%%

[CTCube,CTCubePos,CTThickness]=cubegenerator(PlanCTInfo{1,1},PlanCTLoca);
[CBCTCube,CBCTCubePos,CBCTThickness]=cubegenerator(PreCBCTInfo{1,1},PreCBCTLoca);

subplot(plotrow,plotcol,folderi),
plotcube(CTCube,CTCubePos,0.5,[1 0 0])
plotcube(CBCTCube,CBCTCubePos,0.5,[0 1 0])
title(PlanCTInfo{1, 1}.PatientID);
%%

end
%%
function [CTCube,CTCubePos,CTThickness]=cubegenerator(Sinfo,ScanLoca)
CTCubePos(1)=Sinfo.ImagePositionPatient(1);
CTCubePos(2)=Sinfo.ImagePositionPatient(2);
CTCubePos(3)=Sinfo.ImagePositionPatient(3);
CTCubePos=double(CTCubePos);
CTCube(1)=Sinfo.Rows;
CTCube(2)=Sinfo.Columns;
CTCube(3)=Sinfo.SliceThickness*length(ScanLoca);
CTCube=double(CTCube);
CTThickness=double(Sinfo.SliceThickness);
end
