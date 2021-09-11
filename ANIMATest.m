%% ANIMA Validation Test
clc;clear;close all;
%%
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
addpath('/home/arun/Documents/MATLAB/MaWS');
addpath('/home/arun/Documents/MATLAB/MatOutputs/ANIMAOutput/Run1');
addpath('/home/arun/Documents/MATLAB/ImageDB');
addpath('/home/arun/Documents/MATLAB/ImageOutputs');
addpath('/home/arun/Documents/MATLAB/ImageOutputs/ANIMAOutputs');
addpath('/home/arun/Documents/MATLAB/ImageWorkspace/ANIMAWorkspace');
addpath('/home/arun/src/build/bin');
cd '/home/arun/src/build/bin'
setenv('LD_LIBRARY_PATH','/home/arun/src/build/ITK/lib'); % Mandatory environment for ANIMA
%%
mypath='/home/arun/Documents/MATLAB/ImageDB/Bladder';
patfolder=dir(mypath);
patfolder(1:2,:) = [];
bigfolderlen=length(patfolder)-3;
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
%%
subfolder{1,1}=[];
index=cellfun(@isempty, subfolder) == 0;
subfolder=subfolder(index);
%%
bigfolderlen=length(subfolder);
bigfolderind=1;
% for bigfolderind=1:bigfolderlen
    PlanScanName=subfolder{bigfolderind,1}(1,1).name;
    PathFolder=subfolder{bigfolderind,1}(1,1).folder;
    subfolderlen=length(subfolder{bigfolderind,1});
    directory_name1=char(strcat(PathFolder,'/',PlanScanName));
    
    subfolderind=3;
%     for subfolderind=2:subfolderlen
    TretScanName=subfolder{bigfolderind,1}(subfolderind,1).name;
    directory_name2=char(strcat(PathFolder,'/',TretScanName));
[IP,IT,MP1a,MT1a,Z1a,~,SN1a,SN2a]=loadtwomedimagadatawithmask4woresZ_inbuilt(directory_name1,directory_name2,'CT');
Sname='Structure 2';
TF1x=ismember(SN1a,Sname);
TF2x=ismember(SN2a,Sname);
MP2ac=MP1a(TF1x);
MT2ac=MT1a(TF2x);
MP=MP2ac{1,1};
MT=MT2ac{1,1};


%% sub-routine for ANIMA validation
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
    matworkspacefilename=['/home/arun/Documents/MATLAB/MatOutputs/ANIMAOutput/Run1/',TretScanName,'-run.mat'];
    
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
    %%
    % Dice Score from MATLAB
    NMT=niftiread(filename2Minter);
    IOT=niftiread(transoutputfilename);  
    NIO=niftiread(outputfilename);
    NIO=int16(NIO);
    ITRang=double(max(max(max(IT)))-min(min(min(IT))));
    [resm1,resn1,resp1]=size(NMT);
    CI=zeros(resm1,resn1,resp1);
    MVC=zeros(resp1,1);
    MVT=zeros(resp1,1);
    MVI=zeros(resp1,1);
    MAD=zeros(resp1,1);
    for ind1=1:resp1
    CI(:,:,ind1)=NMT(:,:,ind1)&IOT(:,:,ind1);
    MVI(ind1,1)=length(find(CI(:,:,ind1)));
    MVC(ind1,1)=length(find(IOT(:,:,ind1)));
    MVT(ind1,1)=length(find(NMT(:,:,ind1)));
    
    MAD(ind1,1)=sum(sum(imabsdiff(NIO(:,:,ind1),IT(:,:,ind1))));
    end
    VC=sum(MVC);
    VI=sum(MVI);
    VT=sum(MVT);
    DiceMatlab=(2*VI)/(VT+VC);
    
    MAD1=sum(MAD)/(resm1*resn1*resp1*ITRang);
%     save(matworkspacefilename);
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