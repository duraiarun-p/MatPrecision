clc;clear;close all;
% profile on;
% addpath('/home/s1785969/RDS/MATLAB/DAtoolbox');
% addpath('/home/s1785969/RDS/MATLAB/BlockMatchingAlgoMPEG');
% addpath('/home/s1785969/Documents/MATLAB/MotionDynamics');
addpath('/home/arun/Documents/PyWSPrecision/dblm');
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
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
costV=65537;
m=500;
n=500;
indSize=3;
mbSize=5;
mbSize_N=mbSize*mbSize;
p=5;
ZL=size(IP1a);ZL=ZL(3);

ind=1;
    indzT=max(1-ind,-floor(indSize/2)):min(ZL-ind,floor(indSize/2));
    indL=length(indzT); %Dynamic Z-direction slice choosing
    indz=ind+indzT; % Accessing Z location as index
    Zloc=Z1a(indz);
    ZLocV=(Zloc-Zloc(indzT==0))*-10;

    imgTreatment=IT1a(:,:,ind);
    imgPlanning=IP1a(:,:,indz);

it=imresize(imgTreatment,[m n]);
mp=int16(zeros(m,n,indL));
    for si=1:indL
    mp(:,:,si)=imresize(imgPlanning(:,:,si),[m n]);
    % masTreatment=MT2(:,:,ind);
    % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
    end
imgTreatment=it;
imgPlanning=mp;
clear it mp
 % Converting the Zloc into motion vector and changing the direction from actual Z co-ordinate
%%
% for i = 1 : mbSize : m-(blkfactor*mbSize)+1
%     for j = 1 : mbSize : n-(blkfactor*mbSize)+1
%         currentBlk=imgTreatment(i:i+mbSize-1,j:j+mbSize-1);
%     for zi=1:indL
lambda1A=0.1;
i=1;
j=1;
% zi=1;
prow=(2*p)+1;
costs1 = ones(prow, prow, indL) * costV;
dist = zeros(prow, prow, indL);
refcell=cell(prow, prow, indL);
currentBlk=imgTreatment(i:i+mbSize-1,j:j+mbSize-1);
   for zi=1:indL 
        imgPlanning1=imgPlanning(:,:,zi);
        for m1 = max(1-i,-p) : min(m+1-mbSize-i,p)
            refBlkVer = i + m1;   % m/Vert co-ordinate for ref block            
            imgPsubset1=imgPlanning1(refBlkVer:refBlkVer+mbSize-1,:);
            br = floor((refBlkVer-1)/mbSize)+1;
            for n1 = max(1-j,-p) : min(n+1-mbSize-j,p)
                refBlkHor = j + n1;   % n/Horizontal co-ordinate    
                bc = floor((refBlkHor-1)/mbSize)+1;
                refBlk1=imgPsubset1(((refBlkHor-1)*mbSize+1):(refBlkHor+mbSize-1)*mbSize);
                err1=sum(abs(currentBlk(:)-refBlk1(:)));
                dist(m1+p+1,n1+p+1,zi) = (lambda1A*round(norm([i j]-[refBlkVer refBlkHor])+abs(ZLocV(zi))));
%                 costs1(m1+p+1,n1+p+1,zi) = (err1);
                costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N);
                refcell{m1+p+1,n1+p+1,zi}=refBlk1;
                
%                 costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N)+ dist(m1+p+1,n1+p+1,zi);           
            end
        end
    end 
%    end
% end


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

