%% Visualize DVF as heat map
clc;clear;close all;
addpath('/home/arun/Documents/MATLAB/ImageDB/');
%%
% 
% %%
mbSize=5;m=500;n=500;LenY=round(m/mbSize);LenX=LenY; mbCount=LenX*LenY;
Qz=zeros(mbCount,1);
[Qx,Qy]=quivervector(mbSize,m,n,LenX,LenY);
% 
% a=randi([-5,5],[mbCount,1]);
% b=randi([-5,5],[mbCount,1]);
% c=randi([-1,1],[mbCount,1]);
% % Simulated motion vectors
% mV1=[a,b,c];
% 
% M=sqrt(a.^2+b.^2+c.^2);
% M1=reshape(M,[LenX,LenY]);

mV1=xlsread('motionvector.xlsx');

mc=mV1(:,2);
mcr=reshape(mc,[LenX,LenY]);
% % mC=zeros(LenX,LenY);
% start1=1:100:10000;
% stop1=100:100:10000;
% for i=1:LenX
%     si1=start1(i);
%     si2=stop1(i);
%     mC(:,i)=mc(si1:si2);
% end
Dis=zeros(mbCount,1);
for i=1:mbCount
%     Qy(:,1),Qx(:,1),Qz(:,1),Qy+mV1(:,2),Qx+mV1(:,1),mV1(:,3)/3
    Dis(i,1)=mV1(i,2)^2+mV1(i,2)^2+(mV1(i,3)/3)^2;
end

indSize=max(mV1(:,5));
Zloc=max(mV1(:,3));
DVF=zeros([m,n,indSize]);
mbCount = 1;
    for mi = 1:mbSize:m-mbSize+1
    for mj = 1:mbSize:n-mbSize+1      
        % dy is row(vertical) index
        % dx is col(horizontal) index
        % this means we are scanning in order       
        dy = mV1(mbCount,1);
        dx = mV1(mbCount,2);
        refBlkVerm = mi + dy;
        refBlkHorm = mj + dx;
        refBlkDeptm = mV1(mbCount,5);
%         Dis=sqrt((dy-5)^2+(dx-5)^2+(mV1(mbCount,3)-Zloc)^2);
        Dis=sqrt((dy-5)^2+(dx-5)^2);
%         Dis=(Qy(mbCount)-dy)^2-(Qx(mbCount)-dx)^2+(mV1(mbCount,3)-Zloc)^2;
        DVF(refBlkVerm:refBlkVerm+mbSize-1, refBlkHorm:refBlkHorm+mbSize-1,refBlkDeptm)=Dis;
            mbCount = mbCount + 1;
    end
    end

figure(85);
qplota=quiver3(Qy(:,1),Qx(:,1),Qz(:,1),Qy+mV1(:,2),Qx+mV1(:,1),mV1(:,3)/3);
qplota.Color='b';
qplota.LineWidth=1;
hold on;
% surf(Qy(:,1),Qx(:,1),M1);
hold off;
figure(86);
for ploti=1:indSize
    subplot(3,2,ploti);
    imagesc(DVF(:,:,ploti));
    colormap jet;
    colorbar;
    title(num2str(ploti));
end
%%
%% Motion Compensated Image Generation
function imageComp=motioncomp3d_1_inbuilt(masPlanning,mV1,m,n,indL,mbSize)
% [sr,sc,~]=size(masPlanning);
mp=zeros(m,n,indL);
    for si=1:indL
    mp(:,:,si)=imresize(masPlanning(:,:,si),[m n]);
    % masTreatment=MT2(:,:,ind);
    % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
    end
masPlanning=mp;
imageComp=zeros(m,n);
mbCount = 1;
    for mi = 1:mbSize:m-mbSize+1
    for mj = 1:mbSize:n-mbSize+1      
        % dy is row(vertical) index
        % dx is col(horizontal) index
        % this means we are scanning in order       
        dy = mV1(mbCount,1);
        dx = mV1(mbCount,2);
        refBlkVerm = mi + dy;
        refBlkHorm = mj + dx;
        refBlkDeptm = mV1(mbCount,5);
%         fprintf('mi=%d mj =%d mbCount=%d refBlkVerm = %d, refBlkHorm = %d, dy = %d, dx = %d\n',mi,mj,mbCount,refBlkVerm,refBlkHorm,dy,dx);
        %disp(X);
        imageComp(mi:mi+mbSize-1,mj:mj+mbSize-1) = masPlanning(refBlkVerm:refBlkVerm+mbSize-1, refBlkHorm:refBlkHorm+mbSize-1,refBlkDeptm);
            mbCount = mbCount + 1;
    end
    end
% imageComp=imresize(imageComp,[sr sc]);
end

function [Qx,Qy]=quivervector(mbSize,m,n,LenX,LenY)
if m-n==0
X=1:mbSize:n;X=X';
Y=1:mbSize:n;Y=Y';
else
    X=1:mbSize:n;X=X';
    Y=1:mbSize:m;Y=Y';
end
    CX=floor(X+(mbSize/2));%Grid for Quiver plot
    CY=floor(Y+(mbSize/2));%Grid for Quiver plot

% LenX=length(X);%No. of Blocks
% LenY=length(Y);

vecSize=LenX*LenY;%Total No. of Blocks

Qy=zeros(vecSize,1);
Qx=zeros(vecSize,1);
for qj=1:LenX:vecSize
    qi=(qj+LenX)-1;
    qi1=qi/LenX;
    Qy(qj:qi,1)=CY;
    Qx(qj:qi,1)=CX(qi1);
end
end