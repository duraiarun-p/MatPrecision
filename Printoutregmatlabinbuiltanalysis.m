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

start=2;
stop=start;
% stop=2;
% stop=bigfolderlen;

plotrow=round(sqrt(stop));
if (plotrow^2<stop)
    plotcol=plotrow+1;
else
    plotcol=plotrow;
end
% figure(1),
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
end
%% Registration from MATLAB inuil functions
% helperVolumeRegistration(PreCBCTData,PlanCTData);
fixedHeader=PreCBCTInfo{1,1};
movingHeader=PlanCTInfo{1,1};
fixedVolume=PreCBCTData;
movingVolume=PlanCTData;
%%
centerFixed = round(size(fixedVolume)*0.75);
centerMoving = round(size(movingVolume)*0.75);
%%

%%
[optimizer,metric] = imregconfig('multimodal');
%%
Rfixed  = imref3d(size(fixedVolume),fixedHeader.PixelSpacing(2),fixedHeader.PixelSpacing(1),fixedHeader.SliceThickness);
Rmoving = imref3d(size(movingVolume),movingHeader.PixelSpacing(2),movingHeader.PixelSpacing(1),movingHeader.SliceThickness);
%%
optimizer.InitialRadius = 0.004;
optimizer.MaximumIterations = 200;
% [movingRegisteredVolume, tform] = imregister(movingVolume,Rmoving, fixedVolume,Rfixed, 'rigid', optimizer, metric);
geomtform = imregtform(movingVolume,Rmoving, fixedVolume,Rfixed, 'rigid', optimizer, metric);
movingRegisteredVolume = imwarp(movingVolume,Rmoving,geomtform,'bicubic','OutputView',Rfixed);
%%
maskind1=find(contains(ContourName,'Body','IgnoreCase',true));
% maskind2=find(contains(ContourName,'Bladder','IgnoreCase',true));
maskind2=find(contains(ContourName,'Bladder'));
movingmask1=PlanCTMaskcell{maskind1, 1};
movingRegisteredmask1 = imwarp(movingmask1,Rmoving,geomtform,'bicubic','OutputView',Rfixed);
movingmask2=PlanCTMaskcell{maskind2, 1};
movingRegisteredmask2 = imwarp(movingmask2,Rmoving,geomtform,'bicubic','OutputView',Rfixed);

ctsize=size(movingVolume);
rctsize=size(movingRegisteredVolume);

%%
% figure(2),
% subplot(131),
% imshowpair(movingVolume(:,:,centerMoving(3)), fixedVolume(:,:,centerFixed(3)));title('Before');
% subplot(132),
% imshowpair(movingRegisteredVolume(:,:,centerFixed(3)), fixedVolume(:,:,centerFixed(3)));title('After imregister');
% subplot(133),
% % imshowpair(movingRegisteredVolume2(:,:,47), fixedVolume(:,:,47));title('After imregtform-imwarp');
% imshow(movingRegisteredmask(:,:,centerFixed(3)));
%%
fig = uifigure(1);
ax1=uiaxes(fig);
ax1.Position=[10 200 150 150];
ax1.Title.String = 'Plan CT';
ax2=uiaxes(fig);
ax2.Position=[10 10 150 150];
ax2.Title.String = 'CBCT';
ax3=uiaxes(fig);
ax3.Position=[210 200 150 150];
ax3.Title.String = 'Registered Plan CT';
ax4=uiaxes(fig);
ax4.Position=[210 10 150 150];
ax4.Title.String = 'Overlaid Images';
ax5=uiaxes(fig);
ax5.Position=[410 200 150 150];
ax5.Title.String = 'Registered Plan CT Body Mask';
ax6=uiaxes(fig);
ax6.Position=[410 10 150 150];
ax6.Title.String = 'Registered Plan CT Bladder Mask';

sld = uislider(fig,'Position',[350 200 200 3],'ValueChangedFcn',@(sld,event) updateslide1(sld,fixedVolume,movingRegisteredVolume,movingRegisteredmask1,movingRegisteredmask2,ax2,ax3,ax4,ax5,ax6));
sld.Limits = [1 rctsize(3)];

sld1 = uislider(fig,'Position',[10 200 150 3],'ValueChangedFcn',@(sld1,event) updateslide2(sld1,movingVolume,ax1));
sld1.Limits = [1 ctsize(3)];

imshow(movingVolume(:,:,1),[],'parent',ax1);


imgind=1;
imshow(fixedVolume(:,:,imgind),[],'parent',ax2);

imshow(movingRegisteredVolume(:,:,imgind),[],'parent',ax3);
imshowpair(movingRegisteredVolume(:,:,imgind), fixedVolume(:,:,imgind),'parent',ax4);
imshow(movingRegisteredmask1(:,:,imgind),'parent',ax5);
imshow(movingRegisteredmask2(:,:,imgind),'parent',ax6);

function updateslide1(sld,fixedVolume,movingRegisteredVolume,movingRegisteredmask1,movingRegisteredmask2,ax2,ax3,ax4,ax5,ax6)
imgind=round(sld.Value);
imshow(fixedVolume(:,:,imgind),[],'parent',ax2);
ax2.Title.String = ['CBCT ','Zind=',num2str(imgind)];
imshow(movingRegisteredVolume(:,:,imgind),[],'parent',ax3);
imshowpair(movingRegisteredVolume(:,:,imgind), fixedVolume(:,:,imgind),'parent',ax4);
imshow(movingRegisteredmask1(:,:,imgind),'parent',ax5);
imshow(movingRegisteredmask2(:,:,imgind),'parent',ax6);
end
function updateslide2(sld1,movingVolume,ax1)
imgind1=round(sld1.Value);
imshow(movingVolume(:,:,imgind1),[],'parent',ax1);
ax1.Title.String = ['CT ','Zind=',num2str(imgind1)];
end
% figure(4),
% for zi=1:round(centerFixed(3)/0.75)
%     imshow(movingRegisteredmask(:,:,zi));
%     title(num2str(zi));
%     pause(0.1);
% end