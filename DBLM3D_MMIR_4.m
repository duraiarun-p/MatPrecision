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
% stop=start;
% stop=2;
stop=bigfolderlen;
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
% end
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
%% Non rigid Registration stage
tic;
[ITE,MTE]=dblmnonrigidmaskgen(PlanCTR,PreCBCTData,PreCBCTLoca,PlanMask1);
timeNR=toc;
%%
%Datatypechange
PlanCTData=int16(PlanCTData);
PreCBCTData=int16(PreCBCTData);
PlanCTR=int16(PlanCTR);
ITE=int16(ITE);
PlanMask=int16(PlanMask);
PlanMask1=int16(PlanMask1);
MTE=int16(MTE);

%%
%Filenames
pid=PlanCTInfo{1, 1}.PatientID;
CTsize=size(PlanCTData);
CBCTsize=size(PreCBCTData);

pCTfile=[pid,'-pCT.dcm'];
kVCTfile=[pid,'-preCBCT.dcm'];
pCTRfile=[pid,'-pCTrigid.dcm'];
pCTNRfile=[pid,'-pCTnonrigid.dcm'];
pCTmaskfile=[pid,'-pCTmask.dcm'];
pCTmaskrigidfile=[pid,'-pCTmaskrigid.dcm'];
pCTmasknonrigidfile=[pid,'-pCTmasknonrigid.dcm'];


cd(outpath);
mkdir(pid);
cd([outpath,'/',pid]);
mkdir([outpath,'/',pid,'/','pCT']);
cd([outpath,'/',pid,'/','pCT']);
dicomwrite(reshape(PlanCTData,[CTsize(1) CTsize(2) 1 CTsize(3)]),pCTfile,PlanCTInfo{1,1},'CreateMode','copy');

cd([outpath,'/',pid]);
mkdir([outpath,'/',pid,'/','preCBCT']);
cd([outpath,'/',pid,'/','preCBCT']);
dicomwrite(reshape(PreCBCTData,[CBCTsize(1) CBCTsize(2) 1 CBCTsize(3)]), kVCTfile,PreCBCTInfo{1,1},'CreateMode','copy');

cd([outpath,'/',pid]);
mkdir([outpath,'/',pid,'/','pCTrigid']);
cd([outpath,'/',pid,'/','pCTrigid']);
dicomwrite(reshape(PlanCTR,[CBCTsize(1) CBCTsize(2) 1 CBCTsize(3)]), pCTRfile,PlanCTInfo{1,1},'CreateMode','copy');

cd([outpath,'/',pid]);
mkdir([outpath,'/',pid,'/','pCTnonrigid']);
cd([outpath,'/',pid,'/','pCTnonrigid']);
dicomwrite(reshape(ITE,[CBCTsize(1) CBCTsize(2) 1 CBCTsize(3)]), pCTNRfile,PlanCTInfo{1,1},'CreateMode','copy');

cd([outpath,'/',pid]);
mkdir([outpath,'/',pid,'/','pCTmask']);
cd([outpath,'/',pid,'/','pCTmask']);
dicomwrite(reshape(PlanMask,[CTsize(1) CTsize(2) 1 CTsize(3)]), pCTmaskfile,PlanCTInfo{1,1},'CreateMode','copy');

cd([outpath,'/',pid]);
mkdir([outpath,'/',pid,'/','pCTmaskrigid']);
cd([outpath,'/',pid,'/','pCTmaskrigid']);
dicomwrite(reshape(PlanMask1,[CBCTsize(1) CBCTsize(2) 1 CBCTsize(3)]), pCTmaskrigidfile,PlanCTInfo{1,1},'CreateMode','copy');
% 
cd([outpath,'/',pid]);
mkdir([outpath,'/',pid,'/','pCTmasknonrigid']);
cd([outpath,'/',pid,'/','pCTmasknonrigid']);
dicomwrite(reshape(MTE,[CBCTsize(1) CBCTsize(2) 1 CBCTsize(3)]), pCTmasknonrigidfile,PlanCTInfo{1,1},'CreateMode','copy');


cd(mypath);

%% Visualisation
% rctsize=size(PlanCTR);
fig = uifigure(1);
ax1=uiaxes(fig);
ax1.Position=[10 200 150 150];
ax1.Title.String = 'CBCT';
ax2=uiaxes(fig);
ax2.Position=[10 10 150 150];
ax2.Title.String = 'Rigid Registered Plan CT';
ax3=uiaxes(fig);
ax3.Position=[210 200 150 150];
ax3.Title.String = 'Non Rigid Registered Plan CT';
ax4=uiaxes(fig);
ax4.Position=[210 10 150 150];
ax4.Title.String = 'Mask before registration';
ax5=uiaxes(fig);
ax5.Position=[410 200 150 150];
ax5.Title.String = 'Mask after registration';
ax6=uiaxes(fig);
ax6.Position=[410 10 150 150];
ax6.Title.String = 'Mask overlaid';

sld = uislider(fig,'Position',[10 200 150 3],'ValueChangedFcn',@(sld,event) updateslide1(sld,PreCBCTData,PlanCTR,ITE,PlanMask1,MTE,ax1,ax2,ax3,ax4,ax5,ax6));
sld.Limits = [1 CTsize(3)];

% sld1 = uislider(fig,'Position',[10 200 150 3],'ValueChangedFcn',@(sld1,event) updateslide2(sld1,movingVolume,ax1));
% sld1.Limits = [1 ctsize(3)];
imgind=1;
imshow(PreCBCTData(:,:,imgind),[],'parent',ax1);
imshow(PlanCTR(:,:,imgind),[],'parent',ax2);
imshow(ITE(:,:,imgind),[],'parent',ax3);
imshow(PlanMask1(:,:,imgind),[],'parent',ax4);
imshow(MTE(:,:,imgind),[],'parent',ax5);
imshowpair(PreCBCTData(:,:,imgind), MTE(:,:,imgind),'parent',ax6);
% imshow(pCTwH(:,:,imgind),[],'parent',ax4);

% imshow(PlanCTRH(:,:,imgind),[],'parent',ax5);
% imshow(movingRegisteredmask2(:,:,imgind),'parent',ax6);

end
%% Functions in built

function updateslide1(sld,PreCBCTData,PlanCTR,ITE,PlanMask1,MTE,ax1,ax2,ax3,ax4,ax5,ax6)
imgind=round(sld.Value);
ax1.Title.String = ['CBCT ','Zind=',num2str(imgind)];
imshow(PreCBCTData(:,:,imgind),[],'parent',ax1);
imshow(PlanCTR(:,:,imgind),[],'parent',ax2);
imshow(ITE(:,:,imgind),[],'parent',ax3);
imshow(PlanMask1(:,:,imgind),[],'parent',ax4);
imshow(MTE(:,:,imgind),[],'parent',ax5);
imshowpair(PreCBCTData(:,:,imgind), MTE(:,:,imgind),'parent',ax6);
end