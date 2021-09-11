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
    ind=43;
    indzT=max(1-ind,-floor(indSize/2)):min(ZL-ind,floor(indSize/2));
    indL=length(indzT); %Dynamic Z-direction slice choosing
    indz=ind+indzT; % Accessing Z location as index
    Zloc=Z1(indz);
    ZLocV=(Zloc-Zloc(indzT==0))*-10; % Converting the Zloc into motion vector and changing the direction from actual Z co-ordinate
    imgTreatment=IT1(:,:,ind);
    imgPlanning=IP1(:,:,indz);
    % masTreatment=MT2(:,:,ind);
    masPlanning=MP1(:,:,indz);
it=imresize(imgTreatment,[m n]);
mp=zeros(m,n,indL);
    for si=1:indL
    mp(:,:,si)=imresize(imgPlanning(:,:,si),[m n]);
    % masTreatment=MT2(:,:,ind);
    % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
    end
imgTreatment=it;
imgPlanning=mp;
% costV=65537;%Cost Value
% costV=2^16+1;
prow=2*p + 1;
pcol=2*p + 1;
costs1 = ones(prow, prow, indL) * costV;%Cost matrix
dist = zeros(prow, prow, indL);
mbCount = 1;
mbSize_N=mbSize*mbSize;
mV1 = zeros(round(m*n/mbSize^2),7); % 3D motion vector + Z index
blkfactor=1;
for i = 1 : mbSize : m-(blkfactor*mbSize)+1
    for j = 1 : mbSize : n-(blkfactor*mbSize)+1
        currentBlk=imgTreatment(i:i+mbSize-1,j:j+mbSize-1);
    for zi=1:indL
        imgPlanning1=imgPlanning(:,:,zi);
        for m1 = max(1-i,-p) : min(m+1-mbSize-i,p)
            refBlkVer = i + m1;   % m/Vert co-ordinate for ref block            
            imgPsubset1=imgPlanning1(refBlkVer:refBlkVer+mbSize-1,:);
            for n1 = max(1-j,-p) : min(n+1-mbSize-j,p)
                refBlkHor = j + n1;   % n/Horizontal co-ordinate                                
                refBlk1=imgPsubset1(((refBlkHor-1)*mbSize+1):(refBlkHor+mbSize-1)*mbSize);
                err1=sum(abs(currentBlk(:)-refBlk1(:)));
                dist(m1+p+1,n1+p+1,zi) = (lambda1A*round(norm([i j]-[refBlkVer refBlkHor])+abs(ZLocV(zi))));
%                 costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N);
                costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N)+ dist(m1+p+1,n1+p+1,zi);           
            end
        end
    end        
        [zcostmin,Zminloc]=min(costs1(:));
        [zr,zc,zp]=ind2sub(size(costs1),Zminloc);        
        mV1(mbCount,1) = (zr-p-1);    % row1 co-ordinate for the vector
        mV1(mbCount,2) = (zc-p-1);    % col1 co-ordinate for the vector
        mV1(mbCount,3) = ZLocV(zp);
        mV1(mbCount,4) = indz(zp);
        mV1(mbCount,5) = zp;      
        P1=[m1 n1 ZLocV(zi)];
%         P1=[floor(i+mbSize/2) floor(j+mbSize/2) ZLocV(2)];
        P2=[mV1(mbCount,1) mV1(mbCount,2) mV1(mbCount,3)];
        ThetainSlope=dot(P1,P2)/norm(cross(P1,P2));
        ThetainSlope(isnan(ThetainSlope))=0;
%                 if isnan(ThetainSlope)==1 
%                     ThetainSlope=0; 
%                 end
        mV1(mbCount,6)=atand(ThetainSlope);  
%         mV1(mbCount,7)=costV-zcostmin;
        mV1(mbCount,7)=zcostmin;     % Reliability Score
        mbCount = mbCount + 1;
        costs1 = ones(prow, pcol, indL) * costV;
        dist = zeros(prow, prow, indL);
    end
end
%%
% ZL=length(Z1a);
% [row,col,ZL]=size(IP1a);
% IP1ag=gpuArray(IP1a);
% IP1bg=gpuArray(zeros(row,col,ZL));
% for i=1:ZL
%     IP1bg(:,:,i)=(IP1ag(:,:,i)+IP1ag(:,:,min(i+1,ZL)))/2;
% end
% IP1b=gather(IP1bg);

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