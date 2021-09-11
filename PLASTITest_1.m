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
plastiworkpath='/home/arun/Documents/MATLAB/PlastTestOps/';
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
% Excluding Validation Data blar1
subfolder{1,1}=[];
index=cellfun(@isempty, subfolder) == 0;
subfolder=subfolder(index);
bigfolderlen=length(subfolder);
T=cell(1);
filecount=0;
bigfolderind=1;
% for bigfolderind=1:bigfolderlen % End Loop in the future
    PlanScanName=subfolder{bigfolderind,1}(1,1).name;
    PathFolder=subfolder{bigfolderind,1}(1,1).folder;
    subfolderlen=length(subfolder{bigfolderind,1});
    directory_name1=char(strcat(PathFolder,'/',PlanScanName));    
    subfolderind=2;
%     for subfolderind=2:subfolderlen % End Loop in the future
    TretScanName=subfolder{bigfolderind,1}(subfolderind,1).name;
    directory_name2=char(strcat(PathFolder,'/',TretScanName));

cd(plastiworkpath);
    
[DiceC,DiceCn,MADC,TimeC]=plastidemons(directory_name1,directory_name2,plastiworkpath,PlanScanName,TretScanName);    
%%
function [DiceMatlab,DicePlasti,MAD1,comptime]=plastidemons(directory_name1,directory_name2,plastiworkpath,PlanScanName,TretScanName)
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

MAD1=immse(int16(NIO),IT)/resp1;

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
Close_DicomGui(myhandle1);
Close_DicomGui(myhandle2);
end