clc;clear;close all;
%%
cd '/home/arun/Documents/MATLAB/ImageDB/PrintoutDB';
mypath='/home/arun/Documents/MATLAB/ImageDB/PrintoutDB';
patfolder=dir(mypath);
patfolder(1:2,:) = []; 
bigfolderlen=length(patfolder);

% % bigfolderlen=2;
% for folderi=1:bigfolderlen
% mypath1=[patfolder(folderi).folder,'/',patfolder(folderi).name];
% end
%%
folderi=1;
mypath1=[patfolder(folderi).folder,'/',patfolder(folderi).name];
load(patfolder(folderi).name);
%%
PlanLen=length(PlanCTLoca);
PreCBCTLen=length(PreCBCTLoca);
PosCBCTLen=length(PosCBCTLoca);

PCT=zeros(PlanCTInfo{1, 1}.Rows,PlanCTInfo{1, 1}.Columns,PlanLen);
PrCB=zeros(PreCBCTInfo{1, 1}.Rows,PreCBCTInfo{1, 1}.Columns,PreCBCTLen);
PoCB=zeros(PosCBCTInfo{1, 1}.Rows,PosCBCTInfo{1, 1}.Columns,PosCBCTLen);

for cti=1:PlanLen
    PCT(:,:,cti)=dicomread([PlanCTPath(cti).folder,'/',PlanCTPath(cti).name]);
end

for cbi=1:PreCBCTLen
    PrCB(:,:,cbi)=dicomread([PreCBCTPath(cbi).folder,'/',PreCBCTPath(cbi).name]);
end

for cbi=1:PosCBCTLen
    PoCB(:,:,cbi)=dicomread([PosCBCTPath(cbi).folder,'/',PosCBCTPath(cbi).name]);
end