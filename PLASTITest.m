%% PLASTIMATCH Validation Test
clc;clear;close all;
%%
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
addpath('/home/arun/Documents/MATLAB/MaWS');
addpath('/home/arun/Documents/MATLAB/MAWS_Precision');
addpath('/home/arun/Documents/MATLAB/ImageDB');
addpath('/home/arun/Documents/MATLAB/ImageOutputs/PlastiOutput/');
plastiworkpath='/home/arun/Documents/MATLAB/ImageOutputs/PlastiOutput/';
%% 
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
bigfolderind=3;
% for bigfolderind=1:bigfolderlen
    PlanScanName=subfolder{bigfolderind,1}(1,1).name;
    PathFolder=subfolder{bigfolderind,1}(1,1).folder;
    subfolderlen=length(subfolder{bigfolderind,1});
    directory_name1=char(strcat(PathFolder,'/',PlanScanName));
    
    subfolderind=3;
%     for subfolderind=2:subfolderlen
    TretScanName=subfolder{bigfolderind,1}(subfolderind,1).name;
    directory_name2=char(strcat(PathFolder,'/',TretScanName));
%%
[IP,IT,MP1a,MT1a,Z1,ZL,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
MP2ac=MP1a(TF1x);
MT2ac=MT1a(TF2x);
MP=MP2ac{1,1}; % Planning masks
MT=MT2ac{1,1}; % Treatment masks
clear MP2ac MT2ac
%% Plastimatch command preparation
plastiplanfile=[directory_name1,'/',PlanScanName,'.mhd'];
plastitretfile=[directory_name2,'/',TretScanName,'.mhd'];
plastcommandfile='plasticommand.txt';
commandfilenamepath=[plastiworkpath,plastcommandfile];
plastibspline(plastitretfile,plastiplanfile,plastiworkpath,commandfilenamepath);
%% Plastimatch exectution
cd(plastiworkpath);
syscommand1=['plastimatch ','register ',plastcommandfile];
tic;
system(syscommand1);
comptime=toc;
%% MAD & Dice calculation
NIO=niftiread('warped_img.nii');
NIO=imrotate(NIO,-90);
NIO=fliplr(NIO);
% NIO=imrotate(NIO,-90);
ITr=imrotate(IT,90);
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
    
%     MAD(ind1,1)=sum(sum(imabsdiff(NIO(:,:,ind1),IT(:,:,ind1))));
    end
    VC=sum(MVC);
    VI=sum(MVI);
    VT=sum(MVT);
    DiceMatlab=(2*VI)/(VT+VC);
    
%     MAD1=sum(MAD)/(resm1*resn1*resp1*ITRang);

e=imabsdiff(int16(NIO),IT);
MAD11=sum(sum(sum(e)));

plasttretmaskfile1=[plastiworkpath,'MT.nii'];
niftiwrite(single(NMT),plasttretmaskfile1);

    %Dice score by Plastimatch
syscommand3=['plastimatch ','dice ','MT.nii ','MTE.nii'];
[status, diceinfo]=system(syscommand3);
DicePlasti=diceextract(diceinfo,plastiworkpath);

%%
ind=23;
% fignum=ind;
figure(ind),
subplot(231),
imshow(IT(:,:,ind),[]);title('IT');
subplot(232),
imshow(NIO(:,:,ind),[]);title('NIO');
subplot(233),
imshow(IP(:,:,ind),[]);title('IP');
subplot(234),
imshow(NMT(:,:,ind),[]);title('NMT');
subplot(235),
imshow(IOT(:,:,ind),[]);title('IOT');
subplot(236),
imshow(NMP(:,:,ind),[]);title('NMP');
%%
function plastibspline(fixedfilepath,movingfilepath,plastiworkpath,commandfilenamepath)
% 
fid=fopen(commandfilenamepath,'wt');
fprintf(fid,'[GLOBAL]\n');
fxdimgpath=['fixed=',fixedfilepath,'\n'];
fprintf(fid,fxdimgpath);
movimgpath=['moving=',movingfilepath,'\n'];
fprintf(fid,movimgpath);
outimgpath=['img_out=',plastiworkpath,'warped_img.nii','\n'];
fprintf(fid,outimgpath);
fprintf(fid,['xform_out=bspline_coefficients.txt','\n']);
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