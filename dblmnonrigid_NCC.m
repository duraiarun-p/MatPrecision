function ITE=dblmnonrigid_NCC(PlanCTR,PreCBCTData,PreCBCTLoca)

m=500;
n=500;
indSize=3; % IndSize must always be odd for the symmetrical Z-directional slices
mbSize=5; % Block Size
p=5; % Search Parameter
% lambda1A=0.25;% Distance Regularisation
lambda1A=1;
% lambda3=0.005;
lambda3=0.00001;
costV=2^16+1;
% lambda1A=[0.25 0.5 1 5 10];
% lambda3=[0.25 0.5 1 5 10];
lambda1A=lambda1A';
lambda3=lambda3'; % Orientation Regularisation
% LamL=length(lambda3);
LenX=m/mbSize; % Block col index
LenY=n/mbSize;
%% block matching sub-routine
[~,~,ZL]=size(PlanCTR);

IP1=imresize3(PlanCTR, [m n ZL]);
IT1=imresize3(PreCBCTData, [m n ZL]);
% CI1=false(m,n,ZL);
% MVC1=zeros(ZL,1);
% MVT1=zeros(ZL,1);
% MVI1=zeros(ZL,1);
% MTC1=false(m,n,ZL);
% CI2=false(m,n,ZL);
% MVC2=zeros(ZL,1);
% MVT2=zeros(ZL,1);
% MVI2=zeros(ZL,1);
% MTC2=false(m,n,ZL);
% E=zeros(ZL,1);
% E1=zeros(ZL,1);
% C=zeros(ZL,1);
% ITE=zeros(resm,resn,ZL);
ITE=zeros(m,n,ZL);
tic;
% ind=1;
    parfor ind=1:ZL
%     C(ind,1)=corr2(IP1(:,:,ind),IT1(:,:,ind));   
IP11=IP1;
Z1=PreCBCTLoca;
    indzT=max(1-ind,-floor(indSize/2)):min(ZL-ind,floor(indSize/2));
    indL=length(indzT); %Dynamic Z-direction slice choosing
    indz=ind+indzT; % Accessing Z location as index
    Zloc=Z1(indz);
    ZLocV=(Zloc-Zloc(indzT==0)); % Converting the Zloc into motion vector and changing the direction from actual Z co-ordinate
    imgTreatment=IT1(:,:,ind);
    imgPlanning=IP11(:,:,indz);
    % masTreatment=MT2(:,:,ind);
%     masPlanning=MP2(:,:,indz);
    mV1=BL3DFirstPass_stripped_3_inbuilt(imgTreatment,imgPlanning,mbSize,p,ZLocV,indz,indL,lambda1A,m,n,costV); % 1st Pass
    mV2=BL3DSecndPass_stripped_2_inbuilt(imgTreatment,imgPlanning,mV1,mbSize,p,ZLocV,indz,indL,lambda1A,lambda3,m,n,LenX,LenY,costV); % 2nd Pass
%     MTE1=motioncomp3d_1_inbuilt(masPlanning,mV1,m,n,indL,mbSize);% Mask prediction - 1st Pass
%     MTE2=motioncomp3d_1_inbuilt(masPlanning,mV2,m,n,indL,mbSize);% Mask Prediction - 2nd Pass
%     MTE1=dustthemask(MTE1);% Mask Correction
%     MTE2=dustthemask(MTE2);
    ITEsl=motioncomp3d_1_inbuilt(imgPlanning,mV2,m,n,indL,mbSize);
    % ITE=motioncomp3d_1(imgPlanning,mV2,m,n,indL,mbSize);
%     E(ind,1)=sum(mV1(:,7));
%     E1(ind,1)=sum(mV2(:,7));
%                 MTC1(:,:,ind)=logical(MTE1);
%                 MT2m1=imresize(MT2(:,:,ind),[m n]);              
%                 CI1(:,:,ind)=MT2m1&MTC1(:,:,ind);            
%                 MVI1(ind,1)=length(find(CI1(:,:,ind)));
%                 MVC1(ind,1)=length(find(MTC1(:,:,ind)));
%                 MVT1(ind,1)=length(find(MT2m1));                
%                 MTC2(:,:,ind)=logical(MTE2);
%                 MT2m2=imresize(MT2(:,:,ind),[m n]);              
%                 CI2(:,:,ind)=MT2m2&MTC2(:,:,ind);            
%                 MVI2(ind,1)=length(find(CI2(:,:,ind)));
%                 MVC2(ind,1)=length(find(MTC2(:,:,ind)));
%                 MVT2(ind,1)=length(find(MT2m2));
%                 ITE(:,:,ind)=imresize(ITEsl,[resm resn]);
                    ITE(:,:,ind)=ITEsl;
    end 
end
%     timeNR=toc;
%% Functions in built


%% Block matching - First Pass
function mV1=BL3DFirstPass_stripped_3_inbuilt(imgTreatment,imgPlanning,mbSize,p,ZLocV,indz,indL,lambda1A,m,n,costV)
% [m,n]=size(imgTreatment);
% it=imresize(imgTreatment,[m n]);
% mp=zeros(m,n,indL);
%     for si=1:indL
%     mp(:,:,si)=imresize(imgPlanning(:,:,si),[m n]);
%     % masTreatment=MT2(:,:,ind);
%     % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
%     end
% imgTreatment=it;
% imgPlanning=mp;
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
% for zi=1:indL
%     imgPlanning1=imgPlanning(:,:,zi);
% for i = 1 : mbSize : m-(blkfactor*mbSize)+1
%     for j = 1 : mbSize : n-(blkfactor*mbSize)+1
for j = 1 : mbSize : n-(blkfactor*mbSize)+1
    for i = 1 : mbSize : m-(blkfactor*mbSize)+1
        currentBlk=imgTreatment(i:i+mbSize-1,j:j+mbSize-1);
    for zi=1:indL
%     imgPlanning1=imgPlanning(:,:,zi);
        
        for m1 = max(1-i,-p) : min(m+1-mbSize-i,p)
            refBlkVer = i + m1;   % m/Vert co-ordinate for ref block            
            imgPsubset1=imgPlanning(refBlkVer:refBlkVer+mbSize-1,:,zi);
            for n1 = max(1-j,-p) : min(n+1-mbSize-j,p)
                refBlkHor = j + n1;   % n/Horizontal co-ordinate                                
                refBlk1=imgPsubset1(((refBlkHor-1)*mbSize+1):(refBlkHor+mbSize-1)*mbSize);
%                 err1=sum(abs(currentBlk(:)-refBlk1(:)));
                err1=1-max(max(corrcoef(currentBlk(:),refBlk1(:))));
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
%                 if 
                    ThetainSlope(isnan(ThetainSlope))=0;
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
end


% end
%% Block matching - Second Pass
function mV2=BL3DSecndPass_stripped_2_inbuilt(imgTreatment,imgPlanning,mV1,mbSize,p,ZLocV,indz,indL,lambda1A,lambda3,m,n,LenX,LenY,costV)
% it=imresize(imgTreatment,[m n]);
% mp=zeros(m,n,indL);
%     for si=1:indL
%     mp(:,:,si)=imresize(imgPlanning(:,:,si),[m n]);
%     % masTreatment=MT2(:,:,ind);
%     % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
%     end
% imgTreatment=it;
% imgPlanning=mp;
Vthe=mV1(:,6);
Vthe1=reshape(Vthe,[LenX,LenY]);
% OmegTh=omegdeiff(Vthe1,LenX,LenY);
OmegTh=omegdiif_DrDave(Vthe1,LenX,LenY);
% BlocOd=omegdiif_DrDave(omeg,LenX,LenY)
% costV=65537;%Cost Value
% costV=2^16+1;%Cost Value
prow=2*p + 1;
pcol=2*p + 1;
costs1 = ones(prow, prow, indL) * costV;%Cost matrix
dist = ones(prow, prow, indL);
blkfactor=1;
mbSize_N=mbSize*mbSize;
mV2 = zeros(round(m*n/mbSize^2),7); % 3D motion vector + Z index
mbCount=1;
%     for i = 1 : mbSize : m-(blkfactor*mbSize)+1
%     for j = 1 : mbSize : n-(blkfactor*mbSize)+1
    for j = 1 : mbSize : n-(blkfactor*mbSize)+1
    for i = 1 : mbSize : m-(blkfactor*mbSize)+1
        currentBlk=imgTreatment(i:i+mbSize-1,j:j+mbSize-1);
    for zi=1:indL
%         imgPlanning1=imgPlanning(:,:,zi);
        for m1 = max(1-i,-p) : min(m+1-mbSize-i,p)
            refBlkVer = i + m1;   % m/Vert co-ordinate for ref block 
            br = floor((refBlkVer-1)/mbSize)+1; %row vector for Omeg matrix
            imgPsubset1=imgPlanning(refBlkVer:refBlkVer+mbSize-1,:,zi);
            for n1 = max(1-j,-p) : min(n+1-mbSize-j,p)
                refBlkHor = j + n1;   % n/Horizontal co-ordinate  
                bc = floor((refBlkHor-1)/mbSize)+1;%col vector
                refBlk1=imgPsubset1(((refBlkHor-1)*mbSize+1):(refBlkHor+mbSize-1)*mbSize);
%                 err1=sum(abs(currentBlk(:)-refBlk1(:)));
                err1=1-max(max(corrcoef(currentBlk(:),refBlk1(:))));
                dist(m1+p+1,n1+p+1,zi) = (lambda1A*round(norm([i j]-[refBlkVer refBlkHor])+abs(ZLocV(zi))));
                  costs1(m1+p+1,n1+p+1,zi) = (err1 / mbSize_N) + dist(m1+p+1,n1+p+1,zi) + (lambda3*OmegTh(br,bc));
            end
        end
    end        
        [zcostmin,Zminloc]=min(costs1(:));
        [zr,zc,zp]=ind2sub(size(costs1),Zminloc);     
        mV2(mbCount,1) = (zr-p-1);    % row1 co-ordinate for the vector
        mV2(mbCount,2) = (zc-p-1);    % col1 co-ordinate for the vector
        mV2(mbCount,3) = ZLocV(zp);
        mV2(mbCount,4) = indz(zp);
        mV2(mbCount,5) = zp;       
%         mV2(mbCount,7)=costV-zcostmin;
        mV2(mbCount,7)=zcostmin; % Reliability Score 
        mbCount = mbCount + 1;
        costs1 = ones(prow, pcol, indL) * costV;
        dist = zeros(prow, prow, indL);
    end
    end 
end
%% Motion Compensated Image Generation
function imageComp=motioncomp3d_1_inbuilt(masPlanning,mV1,m,n,~,mbSize)
% [sr,sc,~]=size(masPlanning);
% mp=zeros(m,n,indL);
%     for si=1:indL
%     mp(:,:,si)=imresize(masPlanning(:,:,si),[m n]);
%     % masTreatment=MT2(:,:,ind);
%     % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
%     end
% masPlanning=mp;
imageComp=zeros(m,n);
mbCount = 1;
    for mj = 1:mbSize:n-mbSize+1
    for mi = 1:mbSize:m-mbSize+1      
        % dy is row(vertical) index
        % dx is col(horizontal) index
        % this means we are scanning in order       
        dy = mV1(mbCount,1);
        dx = mV1(mbCount,2);
        refBlkVerm = mi + dy;
        refBlkHorm = mj + dx;
        refBlkDeptm = mV1(mbCount,5);
        imageComp(mi:mi+mbSize-1,mj:mj+mbSize-1) = masPlanning(refBlkVerm:refBlkVerm+mbSize-1, refBlkHorm:refBlkHorm+mbSize-1,refBlkDeptm);
            mbCount = mbCount + 1;
    end
    end
% imageComp=imresize(imageComp,[sr sc]);
end
function BlocOd=omegdiif_DrDave(omeg,LenX,LenY)

BlocOd=zeros(LenX,LenY);
for dispX=-1:1
    for dispY=-1:1
        if (dispX==0) && (dispY==0)
            continue
        end
        BlocOd(max(1-dispY,1):(end-max(dispY,0)), ...
            max(1-dispX,1):(end-max(dispX,0))) = ...
            BlocOd(max(1-dispY,1):(end-max(dispY,0)), ...
            max(1-dispX,1):(end-max(dispX,0))) + ...
            abs(omeg(max(1-dispY,1):(end-max(dispY,0)), ...
            max(1-dispX,1):(end-max(dispX,0))) - ...
            omeg(max(1+dispY,1):(end+min(dispY,0)), ...
            max(1+dispX,1):(end+min(dispX,0))));
    end
end

end

%% Mask Cleaning/Smoothening/Dusting
% function Mn=dustthemask(MTE)
% se = strel('disk',4);
% Mn=imopen(MTE,se);
% Mn=imclose(Mn,se);
% Mn=imfill(Mn,'holes');
% end
