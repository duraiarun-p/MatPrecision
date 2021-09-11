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
%%
maskind1=find(strcmpi(ContourName,'Bladder'));
PlanMask=PlanCTMaskcell{maskind1,1};
%% Rigid Registration stage
tic;
[geomtform,Rfixed,Rmoving,optimizer,metric]=dblmrigidim(PreCBCTData,PlanCTData,PreCBCTInfo{1,1},PlanCTInfo{1,1});
timeR=toc;
PlanCTR = imwarp(PlanCTData,Rmoving,geomtform,'bicubic','OutputView',Rfixed);%Warping scan 
PlanMask1 = imwarp(PlanMask,Rmoving,geomtform,'bicubic','OutputView',Rfixed);%Warping mask
%% Mutual information per block calculation & time Uncomment once for analysis
% start=220;
% stop=224;
% Ip=int16(PlanCTR(start:stop,start:stop,34));
% It=int16(PreCBCTData(start:stop,start:stop,34));
% CBCTsiz=size(Ip);
% NumPixels=CBCTsiz(1)*CBCTsiz(2);
% 
% tic;
% [Hp,dp]=imhist(Ip);
% [Ht,dt]=imhist(It);
% Hp=Hp(Hp~=0);
% % dp=dp(Hp~=0); % Works too
% Ht=Ht(Ht~=0);
% Hp=Hp/sum(Hp);
% Ht=Ht/sum(Ht);
% Sp=-sum(Hp.*log2(Hp));
% St=-sum(Ht.*log2(Ht));
% % indrow = Ip(:)+1; indcol = It(:)+1; %// should same size indrow 
% % indrow = double(Ip(:)) + 1; indcol = double(It(:)) + 1; 
% jointprob = accumarray([Ip(:)+1 It(:)+1], 1); 
% % Ipr=unique(Ip(:));
% % Itr=unique(It(:));
% % jointprob = accumarray([Ipr Itr], 1); 
% jointprob = jointprob / NumPixels; 
% jointprob = jointprob(jointprob~=0);
% Spt = -sum(jointprob.*log2(jointprob));
% E2=Sp+St-Spt;
% timeRMI=toc;
%%
tic;
ITE=dblmnonrigid_SSI_1(PlanCTR,PreCBCTData,PreCBCTLoca);
timeNR=toc;
%%
% m=500;
% n=500;
% indSize=3; % IndSize must always be odd for the symmetrical Z-directional slices
% mbSize=5; % Block Size
% p=5; % Search Parameter
% % lambda1A=0.25;% Distance Regularisation
% lambda1A=0.1;
% % lambda3=0.005;
% lambda3=0.00001;
% costV=2^16+1;
% % lambda1A=[0.25 0.5 1 5 10];
% % lambda3=[0.25 0.5 1 5 10];
% lambda1A=lambda1A';
% lambda3=lambda3'; % Orientation Regularisation
% % LamL=length(lambda3);
% LenX=m/mbSize; % Block col index
% LenY=n/mbSize;
% %% block matching sub-routine
% [~,~,ZL]=size(PlanCTR);
% 
% IP1=imresize3(PlanCTR, [m n ZL]);
% IT1=imresize3(PreCBCTData, [m n ZL]);
% % CI1=false(m,n,ZL);
% % MVC1=zeros(ZL,1);
% % MVT1=zeros(ZL,1);
% % MVI1=zeros(ZL,1);
% % MTC1=false(m,n,ZL);
% % CI2=false(m,n,ZL);
% % MVC2=zeros(ZL,1);
% % MVT2=zeros(ZL,1);
% % MVI2=zeros(ZL,1);
% % MTC2=false(m,n,ZL);
% % E=zeros(ZL,1);
% % E1=zeros(ZL,1);
% % C=zeros(ZL,1);
% % ITE=zeros(resm,resn,ZL);
% ITE=zeros(m,n,ZL);
% tic;
% ind=1;
% %     parfor ind=1:ZL
% %     C(ind,1)=corr2(IP1(:,:,ind),IT1(:,:,ind));   
% IP11=IP1;
% Z1=PreCBCTLoca;
%     indzT=max(1-ind,-floor(indSize/2)):min(ZL-ind,floor(indSize/2));
%     indL=length(indzT); %Dynamic Z-direction slice choosing
%     indz=ind+indzT; % Accessing Z location as index
%     Zloc=Z1(indz);
%     ZLocV=(Zloc-Zloc(indzT==0)); % Converting the Zloc into motion vector and changing the direction from actual Z co-ordinate
%     imgTreatment=IT1(:,:,ind);
%     imgPlanning=IP11(:,:,indz);