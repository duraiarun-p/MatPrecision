clc;clear;close all;
% profile on;
% addpath('/home/s1785969/RDS/MATLAB/DAtoolbox');
% addpath('/home/s1785969/RDS/MATLAB/BlockMatchingAlgoMPEG');
% addpath('/home/s1785969/Documents/MATLAB/MotionDynamics');
addpath('/home/arun/Documents/PyWSPrecision/dblm');
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
addpath('/home/arun/Documents/MATLAB/MAWS_Precision/MIND-SSC');
% /home/arun/Documents/MATLAB/DicomGUI/Scripting


%%
% directory_name1a='/home/arun/Documents/MATLAB/ImageDB/Bladder/blar2/0777394';
% directory_name2a='/home/arun/Documents/MATLAB/ImageDB/Bladder/blar2/E777394';
directory_name1a='/home/arun/Documents/MATLAB/ImageDB/Bladder/blar1/0428705';
directory_name2a='/home/arun/Documents/MATLAB/ImageDB/Bladder/blar1/E428705';
%%
% [IP1a,IT1a,MP1a,MT1a,Z1a,~,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1a,directory_name2a,'CT');
% Sname='Structure 2';
% TF1x=ismember(SN1a,Sname);
% TF2x=ismember(SN2a,Sname);
% MP2ac=MP1a(TF1x);
% MT2ac=MT1a(TF2x);
% MP2a=MP2ac{1,1};
% MT2a=MT2ac{1,1};
[Z1a,IP1a,IT1a,MP2a,MT2a]=blarloadpycall(directory_name1a,directory_name2a);
%%
IP1=IP1a(6:505,6:505,:);
IT1=IT1a(6:505,6:505,:);
MP2=MP2a(6:505,6:505,:);
MT2=MT2a(6:505,6:505,:);
%%
ZL=length(Z1a);
indSize=3;
ind=39;
    indzT=max(1-ind,-floor(indSize/2)):min(ZL-ind,floor(indSize/2));
    indL=length(indzT); %Dynamic Z-direction slice choosing
    indz=ind+indzT; % Accessing Z location as index
    Zloc=Z1a(indz);
    ZLocV=(Zloc-Zloc(indzT==0))*-10;

    imgTreatment=IT1(:,:,ind);
    imgPlanning=IP1(:,:,indz);
%%
% mind=MIND_descriptor(IT1);
r=0;
mT=MIND_descriptor2D(imgTreatment,r);
mP1=MIND_descriptor2D(imgPlanning(:,:,1),r);
mP2=MIND_descriptor2D(imgPlanning(:,:,2),r);
mP3=MIND_descriptor2D(imgPlanning(:,:,3),r);
%%
ms=size(mT);
d1=imabsdiff(imgTreatment,imgPlanning(:,:,1));
Sm=zeros(ms);

for mi=1:ms(3)
    Sm(:,:,mi)=imabsdiff(mT(:,:,mi),mP1(:,:,mi));
end
Sm1=mean(Sm,3);
Sm2=mean(Sm(:));
%%
mT1a=mean(mT,3);
mP1a=mean(mP1,3);
d2=imabsdiff(mT1a,mP1a);
%%
figure(2);subplot(121),imshow(d2,[]);subplot(122),imshow(d1,[]);
% figure(40);
% imshowpair(Sm1,d1);
figure(1),
% ms=size(mT);
for mi=1:ms(3)
    subplot(2,2,mi),
    imshow(mT(:,:,mi),[]);
    title(num2str(mi));
%     pause(0.1);
end
figure(3),
% ms=size(mT);
for mi=1:ms(3)
    subplot(2,2,mi),
    imshow(mP1(:,:,mi),[]);
    title(num2str(mi));
%     pause(0.1);
end
%%
% figure(1),
% subplot(131),
% imshow(mT(:,:,1),[]);
% subplot(132),
% imshow(mP1(:,:,1),[]);
% subplot(133),
% imshow(d1,[]);