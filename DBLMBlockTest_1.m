clc;clear;close all;

%%
load('/home/arun/Documents/MATLAB/ImageDB/BlockTestscan.mat')
%%
costV=65537;
m=500;
n=500;
indSize=3;
mbSize=5;
mbSize_N=mbSize*mbSize;
p=5;
indL=3;
ZLocV=[3;0;-3];
lambda1A=0.1;
i=256;
j=256;
bi=floor((i-1)/mbSize)+1;
bj=floor((j-1)/mbSize)+1;
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
%             br = floor((refBlkVer-1)/mbSize)+1;
            for n1 = max(1-j,-p) : min(n+1-mbSize-j,p)
                refBlkHor = j + n1;   % n/Horizontal co-ordinate    
%                 bc = floor((refBlkHor-1)/mbSize)+1;
                refBlk1=imgPsubset1(((refBlkHor-1)*mbSize+1):(refBlkHor+mbSize-1)*mbSize);
                err1=sum(abs(currentBlk(:)-refBlk1(:)));
                dist(m1+p+1,n1+p+1,zi) = (lambda1A*round(norm([i j]-[refBlkVer refBlkHor])+abs(ZLocV(zi))));
                costs1(m1+p+1,n1+p+1,zi) = (err1);
%                 costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N);
                refcell{m1+p+1,n1+p+1,zi}=refBlk1;
                
%                 costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N)+ dist(m1+p+1,n1+p+1,zi);           
            end
        end
    end 