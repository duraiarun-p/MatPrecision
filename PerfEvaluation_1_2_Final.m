%% Performance Evaluation of BLMIR with ANIMA & Other
% ANIMA is incorporated but Elastix is yet to be done
clc;clear;close all;
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
addpath('/home/arun/Documents/MATLAB/MaWS');
addpath('/home/arun/Documents/MATLAB/MatOutputs/ANIMAOutput/Run1');
addpath('/home/arun/Documents/MATLAB/ImageDB');
addpath('/home/arun/Documents/MATLAB/ImageOutputs');
addpath('/home/arun/Documents/MATLAB/ImageOutputs/ANIMAOutputs');
addpath('/home/arun/Documents/MATLAB/ImageWorkspace/ANIMAWorkspace');
addpath('/home/arun/src/build/bin');

setenv('LD_LIBRARY_PATH','/home/arun/src/build/ITK/lib'); % Mandatory environment for ANIMA

%% Data directory selection
plastiworkpath='/home/arun/Documents/MATLAB/ImageOutputs/PlastiOutput/';
mypath='/home/arun/Documents/MATLAB/ImageDB/Bladder';
patfolder=dir(mypath);
patfolder(1:2,:) = [];
bigfolderlen=length(patfolder);
subfolder=cell(bigfolderlen,1);

for folderi=1:bigfolderlen  
    subfolder{folderi,1}=dir(strcat(mypath,'/',patfolder(folderi).name));
    Names=extractfield(subfolder{folderi,1},'name');
    Names=Names';
    index=find(contains(Names,'.'));
    subfolder{folderi,1}(index)=[];
    Names(index)=[];
end
% We are not Excluding Validation Data blar1 in this run 
% subfolder{1,1}=[];
% index=cellfun(@isempty, subfolder) == 0;
% subfolder=subfolder(index);
% bigfolderlen=length(subfolder);
T=cell(1);
filecount=0;
% bigfolderind=3;
for bigfolderind=1:bigfolderlen % End Loop in the future
    PlanScanName=subfolder{bigfolderind,1}(1,1).name;
    PathFolder=subfolder{bigfolderind,1}(1,1).folder;
    subfolderlen=length(subfolder{bigfolderind,1});
    directory_name1=char(strcat(PathFolder,'/',PlanScanName));    
%     subfolderind=2;
    for subfolderind=2:subfolderlen % End Loop in the future
    TretScanName=subfolder{bigfolderind,1}(subfolderind,1).name;
    directory_name2=char(strcat(PathFolder,'/',TretScanName));
% Registration Method
    cd '/home/arun/src/build/bin';
%     [DiceC,DiceCn,MADC,NMADC,TimeC]=ANIMAValid(directory_name1,directory_name2);

    [CX,NMADUR,MADUR,~,~,DSC2D,DSC1D,MADD,NMADD,TimeD]=BLMIR(directory_name1,directory_name2,3);
    [CX1,NMADUR1,MADUR1,~,~,DSC2D1,DSC1D1,MADD1,NMADD1,TimeD1]=BLMIR(directory_name1,directory_name2,5);
    [CX2,NMADUR2,MADUR2,~,~,DSC2D2,DSC1D2,MADD2,NMADD2,TimeD2]=BLMIR(directory_name1,directory_name2,7);
   
%     [DiceA,DiceAn,MADA,NMADA,TimeA]=plastibspline(directory_name1,directory_name2,plastiworkpath,PlanScanName,TretScanName);
%     [DiceB,DiceBn,MADB,NMADB,TimeB]=plastidemons(directory_name1,directory_name2,plastiworkpath,PlanScanName,TretScanName);
    
    filecount=filecount+1;
    
%     T{filecount,2}=[DiceA,DiceB,DiceC,DSC2D,DSC2D1,DSC2D2,DiceAn,DiceBn,...
%         DiceCn,DSC1D,DSC1D1,DSC1D2,TimeA,TimeB,TimeC,TimeD,TimeD1,TimeD2,MADA,MADB,MADC,MADD,MADD1,MADD2,...
%         NMADA,NMADB,NMADC,NMADD,NMADD1,NMADD2,MADUR,NMADUR,CX];
    T{filecount,2}=[DSC2D,DSC2D1,DSC2D2,DSC1D,DSC1D1,DSC1D2,TimeD,TimeD1,TimeD2,MADD,MADD1,MADD2,...
        NMADD,NMADD1,NMADD2,MADUR,NMADUR,CX];
    T{filecount,1}=TretScanName;
    cd /home/arun/Documents/MATLAB/ImageOutputs/PerfComp/run6;
    save(num2str(filecount));
    end
end

%% Inbuilt Functions
% Plastimatch Demons
function [DiceMatlab,DicePlasti,MAD1,NMAD1,comptime]=plastidemons(directory_name1,directory_name2,plastiworkpath,PlanScanName,TretScanName)
[~,IT,MP1a,MT1a,~,~,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
MP2ac=MP1a(TF1x);
MT2ac=MT1a(TF2x);
MP=MP2ac{1,1}; % Planning masks
MT=MT2ac{1,1}; % Treatment masks
clear MP2ac MT2ac
plastiplanfile=[directory_name1,'/',PlanScanName,'.mhd'];
plastitretfile=[directory_name2,'/',TretScanName,'.mhd'];
plastcommandfile='plasticommand1.txt';
commandfilenamepath=[plastiworkpath,plastcommandfile];
plastidemonscmdgenerate(plastitretfile,plastiplanfile,plastiworkpath,commandfilenamepath);
cd(plastiworkpath);
syscommand1=['plastimatch ','register ',commandfilenamepath];
tic;
system(syscommand1);
comptime=toc;

NIO=niftiread('warped_img1.nii');
NIO=imrotate(NIO,-90);
NIO=fliplr(NIO);
NMP=imrotate(MP,90);
NMT=imrotate(MT,90);

plastplanmaskfile1=[plastiworkpath,'MP1.nii'];
niftiwrite(single(NMP),plastplanmaskfile1);

syscommand2=['plastimatch ','warp ','--input MP1.nii ','--output-img ','MTE1.nii ','--xf ','demons_vf.mha ','--interpolation ','nn'];
system(syscommand2);

IOT=niftiread('MTE1.nii');
% ITRang=double(max(max(max(IT)))-min(min(min(IT))));
    [resm1,resn1,resp1]=size(IT);
    CI=zeros(resm1,resn1,resp1);
    MVC=zeros(resp1,1);
    MVT=zeros(resp1,1);
    MVI=zeros(resp1,1);
%     MAD=zeros(resp1,1);
    for ind1=1:resp1
    CI(:,:,ind1)=NMT(:,:,ind1)&logical(IOT(:,:,ind1));
    MVI(ind1,1)=length(find(CI(:,:,ind1)));
    MVC(ind1,1)=length(find(IOT(:,:,ind1)));
    MVT(ind1,1)=length(find(NMT(:,:,ind1)));    
%     MAD(ind1,1)=sum(sum(imabsdiff(int16(NIO(:,:,ind1)),IT(:,:,ind1))));
    end
    VC=sum(MVC);
    VI=sum(MVI);
    VT=sum(MVT);
    DiceMatlab=(2*VI)/(VT+VC);  
%     MAD1=sum(MAD)/(resm1*resn1*resp1*ITRang);

% e=imabsdiff(int16(NIO),IT);
% MAD1=sum(sum(sum(e)))/(resm1*resn1*resp1);

MAD1=immse(int16(NIO),IT);
NMAD1=MAD1/sqrt(sqrt(sum(int16(NIO(:))))*sqrt(sum(IT(:))));

plasttretmaskfile1=[plastiworkpath,'MT1.nii'];
niftiwrite(single(NMT),plasttretmaskfile1);

    %Dice score by Plastimatch
syscommand3=['plastimatch ','dice ','MT1.nii ','MTE1.nii'];
[~, diceinfo]=system(syscommand3);
DicePlasti=diceextract1(diceinfo,plastiworkpath);
end
function plastidemonscmdgenerate(fixedfilepath,movingfilepath,plastiworkpath,commandfilenamepath)
fid=fopen(commandfilenamepath,'wt');
fprintf(fid,'[GLOBAL]\n');
fxdimgpath=['fixed=',fixedfilepath,'\n'];
fprintf(fid,fxdimgpath);
movimgpath=['moving=',movingfilepath,'\n'];
fprintf(fid,movimgpath);
outimgpath=['img_out=',plastiworkpath,'warped_img1.nii','\n'];
fprintf(fid,outimgpath);
fprintf(fid,['xform_out=demons_vf.mha','\n']);
fprintf(fid,'[STAGE]\n');
fprintf(fid,'xform=vf\n');
% fprintf(fid,'impl=plastimatch\n');
fprintf(fid,'optim=demons\n');
fprintf(fid,'max_its=100\n');
fprintf(fid,'res=4 4 2\n');
fprintf(fid,'demons_std=1\n');
fprintf(fid,'demons_acceleration=5\n');
fprintf(fid,'demons_homogenization=2\n');
fprintf(fid,'demons_filter_width=3 3 5\n');
fprintf(fid,'[STAGE]\n');
fprintf(fid,'res=2 2 1\n');
fprintf(fid,'demons_std=0.8\n');
fprintf(fid,'demons_homogenization=1\n');
fprintf(fid,'demons_filter_width=3 3 5\n');
fprintf(fid,'[STAGE]\n');
fprintf(fid,'res=1 1 1\n');
fprintf(fid,'demons_std=0.6\n');
fprintf(fid,'demons_homogenization=1\n');
fprintf(fid,'demons_filter_width=3 3 5\n');
% fprintf(fid,'[STAGE]\n');
% fprintf(fid,'res=1 1 1\n');
% fprintf(fid,'demons_std=0.4\n');
% fprintf(fid,'demons_homogenization=1\n');
% fprintf(fid,'demons_filter_width=3 3 5\n');
fclose(fid);
end
function DiceP=diceextract1(diceinfo,plastiworkpath)
dicescorefile=[plastiworkpath,'diceinfo1.txt'];
fid=fopen(dicescorefile,'wt');
fprintf(fid,diceinfo);
fclose(fid);
    frm='%s%f';
    fid=fopen(dicescorefile,'r');
    DicePlasti=textscan(fid,frm,'Delimiter',':');
    fclose(fid);
    DiceP=DicePlasti{1,2}(8,1);
end
%Plastimatch b-spline
function [DiceMatlab,DicePlasti,MAD1,NMAD1,comptime]=plastibspline(directory_name1,directory_name2,plastiworkpath,PlanScanName,TretScanName)
[~,IT,MP1a,MT1a,~,~,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
MP2ac=MP1a(TF1x);
MT2ac=MT1a(TF2x);
MP=MP2ac{1,1}; % Planning masks
MT=MT2ac{1,1}; % Treatment masks
clear MP2ac MT2ac
plastiplanfile=[directory_name1,'/',PlanScanName,'.mhd'];
plastitretfile=[directory_name2,'/',TretScanName,'.mhd'];
plastcommandfile='plasticommand.txt';
commandfilenamepath=[plastiworkpath,plastcommandfile];
plastibsplinecmdgenerate(plastitretfile,plastiplanfile,plastiworkpath,commandfilenamepath);
cd(plastiworkpath);
syscommand1=['plastimatch ','register ',commandfilenamepath];
tic;
system(syscommand1);
comptime=toc;

NIO=niftiread('warped_img.nii');
NIO=imrotate(NIO,-90);
NIO=fliplr(NIO);
NMP=imrotate(MP,90);
NMT=imrotate(MT,90);

plastplanmaskfile1=[plastiworkpath,'MP.nii'];
niftiwrite(single(NMP),plastplanmaskfile1);

syscommand2=['plastimatch ','warp ','--input MP.nii ','--output-img ','MTE.nii ','--xf ','bspline_coefficients.txt ','--interpolation ','nn'];
system(syscommand2);

IOT=niftiread('MTE.nii');
% ITRang=double(max(max(max(IT)))-min(min(min(IT))));
    [resm1,resn1,resp1]=size(IT);
    CI=zeros(resm1,resn1,resp1);
    MVC=zeros(resp1,1);
    MVT=zeros(resp1,1);
    MVI=zeros(resp1,1);
%     MAD=zeros(resp1,1);
    for ind1=1:resp1
    CI(:,:,ind1)=NMT(:,:,ind1)&logical(IOT(:,:,ind1));
    MVI(ind1,1)=length(find(CI(:,:,ind1)));
    MVC(ind1,1)=length(find(IOT(:,:,ind1)));
    MVT(ind1,1)=length(find(NMT(:,:,ind1)));    
%     MAD(ind1,1)=sum(sum(imabsdiff(int16(NIO(:,:,ind1)),IT(:,:,ind1))));
    end
    VC=sum(MVC);
    VI=sum(MVI);
    VT=sum(MVT);
    DiceMatlab=(2*VI)/(VT+VC);  
%     MAD1=sum(MAD)/(resm1*resn1*resp1*ITRang);

% e=imabsdiff(int16(NIO),IT);
% MAD1=sum(sum(sum(e)))/(resm1*resn1*resp1);
MAD1=immse(int16(NIO),IT);
NMAD1=MAD1/sqrt(sqrt(sum(int16(NIO(:))))*sqrt(sum(IT(:))));

plasttretmaskfile1=[plastiworkpath,'MT.nii'];
niftiwrite(single(NMT),plasttretmaskfile1);

    %Dice score by Plastimatch
syscommand3=['plastimatch ','dice ','MT.nii ','MTE.nii'];
[~, diceinfo]=system(syscommand3);
DicePlasti=diceextract(diceinfo,plastiworkpath);
end
function plastibsplinecmdgenerate(fixedfilepath,movingfilepath,plastiworkpath,commandfilenamepath)
fid=fopen(commandfilenamepath,'wt');
fprintf(fid,'[GLOBAL]\n');
fxdimgpath=['fixed=',fixedfilepath,'\n'];
fprintf(fid,fxdimgpath);
movimgpath=['moving=',movingfilepath,'\n'];
fprintf(fid,movimgpath);
outimgpath=['img_out=',plastiworkpath,'warped_img.nii','\n'];
fprintf(fid,outimgpath);
fprintf(fid,['xform_out=bspline_coefficients.txt','\n']);
fprintf(fid,['vf_out=bspline_vf.mha','\n']);
fprintf(fid,'[STAGE]\n');
fprintf(fid,'xform=bspline\n');
fprintf(fid,'impl=plastimatch\n');
fprintf(fid,'optim=steepest\n');
fprintf(fid,'max_its=50\n');
fprintf(fid,'regularization_lambda=0.005\n');
fprintf(fid,'grid_spac=10 10 10\n');
fprintf(fid,'res=4 4 2\n');
fprintf(fid,'convergence_tol=1e-8\n');
fprintf(fid,'mi_histogram_bins=200\n');
fprintf(fid,'[STAGE]\n');
fprintf(fid,'max_its=100\n');
fprintf(fid,'grid_spac=5 5 5\n');
fprintf(fid,'res=2 2 1\n');
fprintf(fid,'convergence_tol=1e-12\n');
fprintf(fid,'[STAGE]\n');
fprintf(fid,'max_its=200\n');
fprintf(fid,'grid_spac=1 1 1\n');
fprintf(fid,'res=1 1 1\n');
fprintf(fid,'convergence_tol=1e-16\n');
fclose(fid);
end
function DiceP=diceextract(diceinfo,plastiworkpath)
dicescorefile=[plastiworkpath,'diceinfo.txt'];
fid=fopen(dicescorefile,'wt');
fprintf(fid,diceinfo);
fclose(fid);
    frm='%s%f';
    fid=fopen(dicescorefile,'r');
    DicePlasti=textscan(fid,frm,'Delimiter',':');
    fclose(fid);
    DiceP=DicePlasti{1,2}(8,1);
end
%% ANIMA function
function [DiceMatlab,DiceAnima,MAD1,NMAD1,comptime]=ANIMAValid(directory_name1,directory_name2)
[IP,IT,MP1a,MT1a,~,~,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
MP2ac=MP1a(TF1x);
MT2ac=MT1a(TF2x);
MP=MP2ac{1,1};
MT=MT2ac{1,1};
    % The ANIMA cannot read mhd files hence re-write has nii files
    filename1inter='/home/arun/Documents/MATLAB/ImageWorkspace/ANIMAWorkspace/IP.nii';
    filename2inter='/home/arun/Documents/MATLAB/ImageWorkspace/ANIMAWorkspace/IT.nii';
    filename1Minter='/home/arun/Documents/MATLAB/ImageWorkspace/ANIMAWorkspace/MP.nii';
    filename2Minter='/home/arun/Documents/MATLAB/ImageWorkspace/ANIMAWorkspace/MT.nii';
    
    niftiwrite(IP,filename1inter);
    niftiwrite(IT,filename2inter);
    niftiwrite(single(MP),filename1Minter);
    niftiwrite(single(MT),filename2Minter);
    
    outputfilename='/home/arun/Documents/MATLAB/ImageOutputs/ANIMAOutputs/IPE.nii';    
    transformfilename='/home/arun/Documents/MATLAB/ImageOutputs/ANIMAOutputs/TransVector.nii';
    transformfilenamexml='/home/arun/Documents/MATLAB/ImageOutputs/ANIMAOutputs/TransVectorX.xml';
    transoutputfilename='/home/arun/Documents/MATLAB/ImageOutputs/ANIMAOutputs/MPE.nii';
    diceoutputfilename='/home/arun/Documents/MATLAB/ImageOutputs/ANIMAOutputs/DiceScores.txt';
%     matworkspacefilename=['/home/arun/Documents/MATLAB/MatOutputs/ANIMAOutput/Run1/',TretScanName,'-run.mat'];
    
    syscommand3=['./animaDenseSVFBMRegistration ','-r ',filename2inter,' -m ',filename1inter,...
    ' -o ',outputfilename,' -O ',transformfilename,' -p 4 -l 0'];
    syscommand4=['./animaTransformSerieXmlGenerator -i ',transformfilename,' -o ',transformfilenamexml];
    syscommand5=['./animaApplyTransformSerie -i ',filename1Minter,' -g ',filename2Minter,...
    ' -t ',transformfilenamexml,' -o ',transoutputfilename];
    % syscommand6=['./animaDiceMeasure -o ',diceoutputfilename,' -J -t ',transoutputfilename,' -r ',filename2Minter];
    syscommand6=['./animaDiceMeasure -o ',diceoutputfilename,' -t ',transoutputfilename,' -r ',filename2Minter];
    tic;
    system(syscommand3)
    comptime=toc;
    system(syscommand4);
    system(syscommand5);
    system(syscommand6);
    
    % Dice Score from ANIMA
    frm='%s%f';
    fileID=fopen(diceoutputfilename,'r');
    DiceAnima=textscan(fileID,frm,'Delimiter',':');
    fclose(fileID);
    DiceAnima=DiceAnima{1,2};
    % Dice Score from MATLAB
    NMT=niftiread(filename2Minter);
    IOT=niftiread(transoutputfilename);  
    NIO=niftiread(outputfilename);
    NIO=int16(NIO);
%     ITRang=double(max(max(max(IT)))-min(min(min(IT))));
    [resm1,resn1,resp1]=size(NMT);
    CI=zeros(resm1,resn1,resp1);
    MVC=zeros(resp1,1);
    MVT=zeros(resp1,1);
    MVI=zeros(resp1,1);
%     MAD=zeros(resp1,1);
    for ind1=1:resp1
    CI(:,:,ind1)=NMT(:,:,ind1)&IOT(:,:,ind1);
    MVI(ind1,1)=length(find(CI(:,:,ind1)));
    MVC(ind1,1)=length(find(IOT(:,:,ind1)));
    MVT(ind1,1)=length(find(NMT(:,:,ind1)));    
%     MAD(ind1,1)=sum(sum(imabsdiff(int16(NIO(:,:,ind1)),IT(:,:,ind1))));
    end
    VC=sum(MVC);
    VI=sum(MVI);
    VT=sum(MVT);
    DiceMatlab=(2*VI)/(VT+VC);  
%     MAD1=sum(MAD)/(resm1*resn1*resp1*ITRang);

% e=imabsdiff(int16(NIO),IT);
% MAD1=sum(sum(sum(e)))/(resm1*resn1*resp1);
MAD1=immse(int16(NIO),IT);
NMAD1=MAD1/sqrt(sqrt(sum(int16(NIO(:))))*sqrt(sum(IT(:))));

end
%% Wrapper BLMIR 3D
function [Cr,NMADUR,MADUR,Rel2,Rel1,DSC2,DSC1,MAD1,NMAD1,time]=BLMIR(directory_name1,directory_name2,indSize)
[IP1a,IT1a,MP1a,MT1a,Z1a,~,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
MP2ac=MP1a(TF1x);
MT2ac=MT1a(TF2x);
MP2a=MP2ac{1,1};
MT2a=MT2ac{1,1};
% [IP1a,IT1a,MP2a,MT2a,Z1a,ZLa]=rearrangeslice(IP1a,IT1a,MP2a,MT2a,Z1a); % Common no. of slices selection
%     if flag==0
%     [sr,sc,sp]=size(IP1a);
%     IP1a=randi(2000,[sr,sc,sp],'int16');
%     MPfe2a=false(sr,sc,sp);
%     end
m=500;
n=500;
% indSize=3; % IndSize must always be odd for the symmetrical Z-directional slices
mbSize=5; % Block Size
p=5; % Search Parameter
% lambda1A=0.0905;% Distance Regularisation
% lambda3=0.0045;
lambda1A=0.1;% Distance Regularisation
lambda3=0.0005;
costV=2^16+1;
% lambda1A=[0.25 0.5 1 5 10];
% lambda3=[0.25 0.5 1 5 10];
% lambda1A=lambda1A';
% lambda3=lambda3'; % Orientation Regularisation
% LamL=length(lambda3);
LenX=m/mbSize; % Block col index
LenY=n/mbSize; % Block row index
[Cr,NMADUR,MADUR,Rel2,Rel1,DSC2,DSC1,MAD1,NMAD1,time]=BLMIR3DDiceFx_24_inbuilt(Z1a,IP1a,IT1a,MP2a,MT2a,indSize,mbSize,p,lambda1A,lambda3,m,n,LenX,LenY,costV);%1 lambda
end

%% Dice score, Reliability score & Block matching for all slices using parfor
function [Cr,NMADUR,MADUR,ECV1,ECV,DSC2,DSC1,MAD1,NMAD1,time]=BLMIR3DDiceFx_24_inbuilt(Z1,IP1,IT1,MP2,MT2,indSize,mbSize,p,lambda1A,lambda3,m,n,LenX,LenY,costV)
ZL=length(Z1); % Total length of slices per scan
CI1=false(m,n,ZL);
MVC1=zeros(ZL,1);
MVT1=zeros(ZL,1);
MVI1=zeros(ZL,1);
MTC1=false(m,n,ZL);
CI2=false(m,n,ZL);
MVC2=zeros(ZL,1);
MVT2=zeros(ZL,1);
MVI2=zeros(ZL,1);
MTC2=false(m,n,ZL);
E=zeros(ZL,1);
E1=zeros(ZL,1);
C=zeros(ZL,1);
% MAD=zeros(ZL,1);
[resm1,resn1,resp1]=size(IT1);
ITEV=zeros(resm1,resn1,resp1);
tic;
    parfor ind=1:ZL
    C(ind,1)=corr2(IP1(:,:,ind),IT1(:,:,ind));   
    indzT=max(1-ind,-floor(indSize/2)):min(ZL-ind,floor(indSize/2));
    indL=length(indzT); %Dynamic Z-direction slice choosing
    indz=ind+indzT; % Accessing Z location as index
    Zloc=Z1(indz);
    ZLocV=(Zloc-Zloc(indzT==0))*-10; % Converting the Zloc into motion vector and changing the direction from actual Z co-ordinate
    imgTreatment=IT1(:,:,ind);
    imgPlanning=IP1(:,:,indz);
    masPlanning=MP2(:,:,indz);
    mV1=BL3DFirstPass_stripped_3_inbuilt(imgTreatment,imgPlanning,mbSize,p,ZLocV,indz,indL,lambda1A,m,n,costV); % 1st Pass
    mV2=BL3DSecndPass_stripped_2_inbuilt(imgTreatment,imgPlanning,mV1,mbSize,p,ZLocV,indz,indL,lambda1A,lambda3,m,n,LenX,LenY,costV); % 2nd Pass
    MTE1=motioncomp3d_1_inbuilt(masPlanning,mV1,m,n,indL,mbSize);% Mask prediction - 1st Pass
    MTE2=motioncomp3d_1_inbuilt(masPlanning,mV2,m,n,indL,mbSize);% Mask Prediction - 2nd Pass
    ITE=motioncomp3d_1_inbuilt(imgPlanning,mV2,m,n,indL,mbSize);% Image Slice Prediction - 2nd Pass
    MTE1=dustthemask(MTE1);% Mask Correction
    MTE2=dustthemask(MTE2);
    E(ind,1)=sum(mV1(:,7));
    E1(ind,1)=sum(mV2(:,7));
                MTC1(:,:,ind)=logical(MTE1);
                MT2m1=imresize(MT2(:,:,ind),[m n]);              
                CI1(:,:,ind)=MT2m1&MTC1(:,:,ind);            
                MVI1(ind,1)=length(find(CI1(:,:,ind)));
                MVC1(ind,1)=length(find(MTC1(:,:,ind)));
                MVT1(ind,1)=length(find(MT2m1));                
                MTC2(:,:,ind)=logical(MTE2);
                MT2m2=imresize(MT2(:,:,ind),[m n]);              
                CI2(:,:,ind)=MT2m2&MTC2(:,:,ind);            
                MVI2(ind,1)=length(find(CI2(:,:,ind)));
                MVC2(ind,1)=length(find(MTC2(:,:,ind)));
                MVT2(ind,1)=length(find(MT2m2));
                ITEV(:,:,ind)=imresize(ITE,[resm1 resn1]);
%                 MAD(ind,1)=sum(sum(imabsdiff(int16(ITE),imresize(imgTreatment,[m n]))));
    end
time=toc;
        VC1=sum(MVC1);
        VT1=sum(MVT1);
        VI1=sum(MVI1);
        DSC1=(2*VI1)/(VT1+VC1);      
        VC2=sum(MVC2);
        VT2=sum(MVT2);
        VI2=sum(MVI2);
        DSC2=(2*VI2)/(VT2+VC2);       
        ECV1=sum(E1)/(ZL*LenX*LenY*costV);
        ECV=sum(E)/(ZL*LenX*LenY*costV);
        Cr=mean(C);
%         ITRang=double(max(max(max(IT1)))-min(min(min(IT1))));
%         MAD1=sum(MAD)/(ZL*m*n*ITRang);

%         e=imabsdiff(int16(ITEV),IT1);
% MAD1=sum(sum(sum(e)))/(resm1*resn1*resp1);
% 
%         E=imabsdiff(IP1,IT1);
% MADUR=sum(sum(sum(E)))/(resm1*resn1*resp1);

MADUR=immse(IP1,IT1);
NMADUR=MADUR/sqrt(sqrt(sum(IP1(:)))*sqrt(sum(IT1(:))));

MAD1=immse(int16(ITEV),IT1);
NMAD1=MAD1/sqrt(sqrt(sum(int16(ITEV(:))))*sqrt(sum(IT1(:))));

end
%% Block matching - First Pass
function mV1=BL3DFirstPass_stripped_3_inbuilt(imgTreatment,imgPlanning,mbSize,p,ZLocV,indz,indL,lambda1A,m,n,costV)
% [m,n]=size(imgTreatment);
it=imresize(imgTreatment,[m n]);
mp=int16(zeros(m,n,indL));
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


end
%% Block matching - Second Pass
function mV2=BL3DSecndPass_stripped_2_inbuilt(imgTreatment,imgPlanning,mV1,mbSize,p,ZLocV,indz,indL,lambda1A,lambda3,m,n,LenX,LenY,costV)
it=imresize(imgTreatment,[m n]);
mp=int16(zeros(m,n,indL));
    for si=1:indL
    mp(:,:,si)=imresize(imgPlanning(:,:,si),[m n]);
    % masTreatment=MT2(:,:,ind);
    % m(:,:,si)=(imresize(masPlanning(:,:,si),[m n]));
    end
imgTreatment=it;
imgPlanning=mp;
Vthe=mV1(:,6);
Vthe1=reshape(Vthe,[LenX,LenY]);
% OmegTh=omegdeiff(Vthe1,LenX,LenY);
OmegTh=omegdiif_DrDave_inbuilt(Vthe1,LenX,LenY);
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
    for i = 1 : mbSize : m-(blkfactor*mbSize)+1
    for j = 1 : mbSize : n-(blkfactor*mbSize)+1
        currentBlk=imgTreatment(i:i+mbSize-1,j:j+mbSize-1);
    for zi=1:indL
        imgPlanning1=imgPlanning(:,:,zi);
        for m1 = max(1-i,-p) : min(m+1-mbSize-i,p)
            refBlkVer = i + m1;   % m/Vert co-ordinate for ref block 
            br = floor((refBlkVer-1)/mbSize)+1; %row vector for Omeg matrix
            imgPsubset1=imgPlanning1(refBlkVer:refBlkVer+mbSize-1,:);
            for n1 = max(1-j,-p) : min(n+1-mbSize-j,p)
                refBlkHor = j + n1;   % n/Horizontal co-ordinate  
                bc = floor((refBlkHor-1)/mbSize)+1;%col vector
                refBlk1=imgPsubset1(((refBlkHor-1)*mbSize+1):(refBlkHor+mbSize-1)*mbSize);
                err1=sum(abs(currentBlk(:)-refBlk1(:)));
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
%% Neighbourhood Summation
function BlocOd=omegdiif_DrDave_inbuilt(omeg,LenX,LenY)
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
        imageComp(mi:mi+mbSize-1,mj:mj+mbSize-1) = masPlanning(refBlkVerm:refBlkVerm+mbSize-1, refBlkHorm:refBlkHorm+mbSize-1,refBlkDeptm);
            mbCount = mbCount + 1;
    end
    end
% imageComp=imresize(imageComp,[sr sc]);
end

%% Mask Cleaning/Smoothening/Dusting
function Mn=dustthemask(MTE)
se = strel('disk',4);
Mn=imopen(MTE,se);
Mn=imclose(Mn,se);
Mn=imfill(Mn,'holes');
end
%% Data Loading DicomGui Scripts
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