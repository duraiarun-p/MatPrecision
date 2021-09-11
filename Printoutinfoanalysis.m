%% PRINTOUT DATA EXTRACTION - FINAL
clc;clear;close all;
% profile on;
% addpath('/home/arun/Documents/MATLAB/DAtoolbox');
% addpath('/home/arun/Documents/MATLAB/BlockMatchingAlgoMPEG');
% addpath('/home/s1785969/Documents/MATLAB/MotionDynamics');
% addpath('/home/arun/Documents/MATLAB/DA_Image_REG');
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');
addpath('/home/arun/Documents/MATLAB/ImageDB/PRINTOUT_TRIAL/Anonymised_Data');
%%

mypath='/home/arun/Documents/MATLAB/ImageDB/PRINTOUT_TRIAL/Anonymised_Data/';
patfolder=dir(mypath);
patfolder(1:2,:) = []; 
bigfolderlen=length(patfolder);
% bigfolderlen=1;
subfolder=cell(bigfolderlen,1);
DataInfocell=cell(1,1);
bigfolderlen=2;
for folderi=1:bigfolderlen
mypath1=[patfolder(folderi).folder,'/',patfolder(folderi).name];
% DataInfocell{folderi,1}=loadprintoudatainbuilt(mypath1);
[PlanCTPath,PlanCTInfo,PlanCTLoca,PreCBCTInfo,PreCBCTPath,PreCBCTLoca,PosCBCTInfo,PosCBCTPath,PosCBCTLoca]=extractprintoutdata_ib(mypath1,'CT');
%%

fileID=randi(length(PreCBCTLoca));
SOPC=PlanCTInfo{fileID, 1}.SOPClassUID;
SOPI=PlanCTInfo{fileID, 1}.SOPInstanceUID;
StudI=PlanCTInfo{fileID, 1}.StudyInstanceUID;
SerieI=PlanCTInfo{fileID, 1}.SeriesInstanceUID;
TranI=PlanCTInfo{fileID, 1}.TransferSyntaxUID;

SOPC1=PreCBCTInfo{fileID, 1}.SOPClassUID;
SOPI1=PreCBCTInfo{fileID, 1}.SOPInstanceUID;
StudI1=PreCBCTInfo{fileID, 1}.StudyInstanceUID;
SerieI1=PreCBCTInfo{fileID, 1}.SeriesInstanceUID;
TranI1=PreCBCTInfo{fileID, 1}.TransferSyntaxUID;

SOPC2=PosCBCTInfo{fileID, 1}.SOPClassUID;
SOPI2=PosCBCTInfo{fileID, 1}.SOPInstanceUID;
StudI2=PosCBCTInfo{fileID, 1}.StudyInstanceUID;
SerieI2=PosCBCTInfo{fileID, 1}.SeriesInstanceUID;
TranI2=PosCBCTInfo{fileID, 1}.TransferSyntaxUID;
%%
cd(mypath1);
RSFiles=dir('RS*.dcm');
filecount=length(RSFiles);
RSinfo=cell(filecount,1);
RSstud=cell(filecount,1);
RSTran=cell(filecount,1);
RSSTim=cell(filecount,2);
    for fi=1:filecount
    RSinfo{fi,1}=dicominfo([RSFiles(fi).folder,'/',RSFiles(fi).name]);
    RSstud{fi,1}=RSinfo{fi,1}.StructureSetLabel;
    RSTran{fi,1}=RSinfo{fi,1}.TransferSyntaxUID;
    RSSTim{fi,1}=RSinfo{fi,1}.StructureSetTime;
    RSSTim{fi,2}=RSinfo{fi,1}.StructureSetDate;
    end
DataInfocell{folderi,1}=RSstud;

end

cd /home/arun/Documents/MATLAB/ImageOutputs;
DataT=cell2table(DataInfocell{1,1});
for folderi=2:bigfolderlen
    Data=cell2table(DataInfocell{folderi,1});
    DataT=[DataT;Data];
end
%%
function [PlanCTPath,PlanCTInfo,PlanCTLoca,PreCBCTInfo,PreCBCTPath,PreCBCTLoca,PosCBCTInfo,PosCBCTPath,PosCBCTLoca]=extractprintoutdata_ib(mypath,Modality)
[CTInfoCell,CTPathCell,CBCTInfocell,CBCTPathcell]=extractprintoutdatainfo_ib(mypath,Modality);

CBCTLen=length(CBCTPathcell);
Preind=CBCTLen-1;
Posind=CBCTLen;
PreCBCTInfo=CBCTInfocell{Preind,1};
PreCBCTPath=CBCTPathcell{Preind,1};
PreCBCTLoca=CBCTInfocell{Preind,2};
PosCBCTInfo=CBCTInfocell{Posind,1};
PosCBCTPath=CBCTPathcell{Posind,1};
PosCBCTLoca=CBCTInfocell{Posind,2};
PlanCTInfo=CTInfoCell{1,1};
PlanCTPath=CTPathCell{1,1};
PlanCTLoca=CTInfoCell{1,2};
end

function [CTInfoCell,CTPathCell,CBCTInfocell,CBCTPathcell]=extractprintoutdatainfo_ib(mypath,Modalityselect)
[CTInfoCell,CTPathCell,CBCSeriesFilesInfoCell,CBCSeriesFilesPathCell]=loadprintoutdata_ib(mypath,Modalityselect);

CBCTLen=length(CBCSeriesFilesInfoCell);
count=1;
for cbct=1:CBCTLen
    ctsublen=length(CBCSeriesFilesInfoCell{cbct,1});
    for cbctsub=1:ctsublen
        if iscell(CBCSeriesFilesInfoCell{cbct,1}{1,cbctsub})==1
            SNo=CBCSeriesFilesInfoCell{cbct,1}{1,cbctsub}{1,1}.SeriesNumber;
            NoS=length(CBCSeriesFilesInfoCell{cbct,1}{1,cbctsub});
            Man=CBCSeriesFilesInfoCell{cbct,1}{1,cbctsub}{1,1}.Manufacturer;
            Tim(count,1)=round(str2double(CBCSeriesFilesInfoCell{cbct,1}{1,cbctsub}{1,1}.SeriesTime)/100);           
            CBCTInfocell{count,1}=CBCSeriesFilesInfoCell{cbct,1}{1,cbctsub};
            CBCTInfocell{count,2}=CBCSeriesFilesInfoCell{cbct,1}{1,cbctsub+2};
            CBCTPathcell{count,1}=CBCSeriesFilesPathCell{cbct,1}{1,cbctsub};
            count=count+1;
        end
    end
end

CBCTLen1=length(CBCTInfocell);
CBCTTime=zeros(CBCTLen1,1);
for cbct=1:CBCTLen1
            Tim=round(str2double(CBCTInfocell{cbct,1}{1,1}.SeriesTime)/10);
            CBCTTime(cbct)=Tim;
end
[CBCTTimeSorted,CBCTTimeIndex]=sort(CBCTTime,'ascend');

CBCTInfocell=CBCTInfocell(CBCTTimeIndex,:);
CBCTPathcell=CBCTPathcell(CBCTTimeIndex);
end

function [CTSeriesFilesInfoCell,CTSeriesFilesPathCell,CBCSeriesFilesInfoCell_1,CBCSeriesFilesPathCell_1]=loadprintoutdata_ib(mypath,ModalitySelect)
patfolder=dir(mypath);
patfolder(1:2)=[];
FoldLen=length(patfolder);
% FoldLen=round(length(patfolder)*.015);
% cd(mypath);
% PosV=zeros(FoldLen,3);
% CBCTIdent=zeros(FoldLen,1);
% SLoc=zeros(FoldLen,1);

Modal=cell(FoldLen,1);
infocell=cell(FoldLen,1);

% imgcell=cell(FoldLen,3);
% parpoolobj=parpool('LocalProfile2',10);
parfor filenumber=1:FoldLen
filename=[patfolder(filenumber).folder,'/',patfolder(filenumber).name];
infocell{filenumber,1}=dicominfo(filename);
Modal{filenumber, 1}=infocell{filenumber,1}.Modality;
end
% delete(parpoolobj);
%% Filtering based on Modality
%Finding all modality information
% [ModalityInFiles,ModalFileStartIndex,~]=unique(Modal);
% ModalFileStopIndex=[ModalFileStartIndex(2:end);FoldLen]-1;
% %Choosing particular modality using string flag
% ModalitySelect='CT';
% ModalitySelectIndex=find(strcmp(ModalityInFiles,ModalitySelect));
ModalityIndex=zeros(FoldLen,1);
for filenumber=1:FoldLen
ModalityIndex(filenumber)=isequal(Modal{filenumber,1},ModalitySelect);
end
%Selecting Dicom files based modality index
% SelectedDicomFiles_Modality=infocell(ModalFileStartIndex(ModalitySelectIndex):ModalFileStopIndex(ModalitySelectIndex));
% SelectedDicomFilesPath=patfolder(ModalFileStartIndex(ModalitySelectIndex):ModalFileStopIndex(ModalitySelectIndex));
SelectedDicomFiles_Modality=infocell(ModalityIndex>0);
SelectedDicomFilesPath=patfolder(ModalityIndex>0);
%% Filtering based on Manufacturer to classifiy CT from CBCT
Manufacture_ID=zeros(length(SelectedDicomFiles_Modality),1);
for filenumber=1:length(SelectedDicomFiles_Modality)
%     if strcmp(infocell{filenumber,1}.Manufacturer,'Philips')==1
    if isfield(infocell{filenumber,1},'RescaleType')==1
        Manufacture_ID(filenumber,:)=0;
    else
        Manufacture_ID(filenumber,:)=1;
    end
end
CTDicomFiles=SelectedDicomFiles_Modality(Manufacture_ID>0);
CTDicomFilespath=SelectedDicomFilesPath(Manufacture_ID>0);
CBCTDicomFiles=SelectedDicomFiles_Modality(Manufacture_ID==0);
CBCTDicomFilespath=SelectedDicomFilesPath(Manufacture_ID==0);

%CBCT Series Number Extraction
CBCTSerie=zeros(length(CBCTDicomFiles),1);
for filenumber=1:length(CBCTDicomFiles)
    CBCTSerie(filenumber,:)=CBCTDicomFiles{filenumber,1}.SeriesNumber;
end
CTSerie=zeros(length(CTDicomFiles),1);
for filenumber=1:length(CTDicomFiles)
    CTSerie(filenumber,:)=CTDicomFiles{filenumber,1}.SeriesNumber;
end

% CBCT Filtering based on Series Number
[CBCTSeriesinFiles,~,~]=unique(CBCTSerie);
CBCTSerLen=length(CBCTSeriesinFiles);
CBCTSerInd=zeros(length(CBCTDicomFiles),length(CBCTSeriesinFiles));
CBCTSersubLen=zeros(1,CBCTSerLen);
CBCSeriesFilesInfoCell=cell(CBCTSerLen,1);
CBCSeriesFilesPathCell=cell(CBCTSerLen,1);
for ctseri=1:CBCTSerLen
      for filenumber=1:length(CBCTDicomFiles)
          if CBCTSerie(filenumber)==CBCTSeriesinFiles(ctseri)
            CBCTSerInd(filenumber,ctseri)=1;
          else
              CBCTSerInd(filenumber,ctseri)=0;
          end
      end
   CBCTSersubLen(1,ctseri)=length(find(CBCTSerInd(:,ctseri)));
   CBCSeriesFilesInfoCell{ctseri,1}=CBCTDicomFiles(CBCTSerInd(:,ctseri)>0);
   CBCSeriesFilesPathCell{ctseri,1}=CBCTDicomFilespath(CBCTSerInd(:,ctseri)>0);
end
% CBCT Filtering based on Time of Acquisition (Series Time Dicom attribute)
CBCSeriesFilesInfoCell_1=cell(CBCTSerLen,1);
CBCSeriesFilesPathCell_1=cell(CBCTSerLen,1);
for ctseri=1:CBCTSerLen
    CBCSeriesTime=zeros(CBCTSersubLen(1,ctseri),1);
    for filenumber=1:CBCTSersubLen(1,ctseri)
        CBCSeriesTime(filenumber)=round(str2double(CBCSeriesFilesInfoCell{ctseri,1}{filenumber,1}.SeriesTime)/100);
    end
    [CBCTSeriesTimeinFiles,~,~]=unique(CBCSeriesTime);
    if length(CBCTSeriesTimeinFiles)>=2
    CBCTSeriesTimeinFiles_select=[CBCTSeriesTimeinFiles(end-1);CBCTSeriesTimeinFiles(end)];
    else
    CBCTSeriesTimeinFiles_select=CBCTSeriesTimeinFiles(end);
    end
    CBCTSeriesTimeInd=zeros(CBCTSersubLen(1,ctseri),1);
        for filenumber=1:CBCTSersubLen(1,ctseri)
            %removing minutes and seconds discrepencies
            if round(str2double(CBCSeriesFilesInfoCell{ctseri,1}{filenumber,1}.SeriesTime)/100)==CBCTSeriesTimeinFiles_select(1)
                CBCTSeriesTimeInd(filenumber)=1;
            elseif round(str2double(CBCSeriesFilesInfoCell{ctseri,1}{filenumber,1}.SeriesTime)/100)==CBCTSeriesTimeinFiles_select(2)
                CBCTSeriesTimeInd(filenumber)=2;
            end       
        end
        CBCTSeriesTimeIndVal=[1;2];
        CBCTemp=CBCSeriesFilesInfoCell{ctseri,1};
        CBCTemp1=CBCSeriesFilesPathCell{ctseri,1};
        CBCSeriesPrePostCell=cell(1,length(CBCTSeriesTimeinFiles_select));
        CBCSeriesPrePostCell_1=cell(1,length(CBCTSeriesTimeinFiles_select));
        for cbctserseli=1:length(CBCTSeriesTimeinFiles_select)
        CBCSeriesPrePostCell{1,cbctserseli}=CBCTemp(CBCTSeriesTimeInd(:,1)==CBCTSeriesTimeIndVal(cbctserseli));
        CBCSeriesPrePostCell_1{1,cbctserseli}=CBCTemp1(CBCTSeriesTimeInd(:,1)==CBCTSeriesTimeIndVal(cbctserseli));
        end
        
        CBCSeriesFilesInfoCell_1{ctseri,1}=CBCSeriesPrePostCell;
        CBCSeriesFilesPathCell_1{ctseri,1}=CBCSeriesPrePostCell_1;
end
% CBCT Filtering & Sorting based on Slice Location
for ctseri=1:CBCTSerLen
    CBCTPrePoLen=length(CBCSeriesFilesInfoCell_1{ctseri,1});
    for cbctprepoi=1:CBCTPrePoLen
        CBCTLoc=zeros(length(CBCSeriesFilesInfoCell_1{ctseri,1}{1,cbctprepoi}),1);
        for filenumber=1:length(CBCSeriesFilesInfoCell_1{ctseri,1}{1,cbctprepoi})
        CBCTLoc(filenumber)=CBCSeriesFilesInfoCell_1{ctseri,1}{1,cbctprepoi}{filenumber,1}.ImagePositionPatient(3);
        CBCSeriesFilesInfoCell_1{ctseri,1}{1,cbctprepoi}{filenumber,2}=CBCTLoc(filenumber);
        end
        [CBCTLoc_sorted,CBCTLoc_sorted_ind]=sort(CBCTLoc,'descend');
        CBCSeriesFilesInfoCell_1{ctseri,1}{1,cbctprepoi}=CBCSeriesFilesInfoCell_1{ctseri,1}{1,cbctprepoi}(CBCTLoc_sorted_ind);
        CBCSeriesFilesPathCell_1{ctseri,1}{1,cbctprepoi}=CBCSeriesFilesPathCell_1{ctseri,1}{1,cbctprepoi}(CBCTLoc_sorted_ind);
        CBCSeriesFilesInfoCell_1{ctseri,1}{1,cbctprepoi+2}=CBCTLoc_sorted;
    end
end


% CT Series Extraction& Filtering
[CTSeriesinFiles,~,~]=unique(CTSerie);
CTSerLen=length(CTSeriesinFiles);
CTSerInd=zeros(length(CTDicomFiles),length(CTSeriesinFiles));
CTSersubLen=zeros(1,CTSerLen);
CTSeriesFilesInfoCell=cell(CTSerLen,2);
CTSeriesFilesPathCell=cell(CTSerLen,1);
for ctseri=1:CTSerLen
      for filenumber=1:length(CTDicomFiles)
          if CTSerie(filenumber)==CTSeriesinFiles(ctseri)
            CTSerInd(filenumber,ctseri)=1;
          else
              CTSerInd(filenumber,ctseri)=0;
          end
      end
   CTSersubLen(1,ctseri)=length(find(CTSerInd(:,ctseri)));
   CTSeriesFilesInfoCell{ctseri,1}=CTDicomFiles(CTSerInd(:,ctseri)>0);
   CTSeriesFilesPathCell{ctseri,1}=CTDicomFilespath(CTSerInd(:,ctseri)>0);
end
% CTSerLen=length(CTSeriesinFiles);
CTLen=zeros(CTSerLen,1);
for ctseri=1:CTSerLen
    CTLen(ctseri)=length(CTSeriesFilesInfoCell{ctseri,1});
end
CTSeriesFilesPathCell(CTLen<=1,:)=[];
CTSeriesFilesInfoCell(CTLen<=1,:)=[];
CTSerLen=length(CTSeriesFilesPathCell);
%% Sorting CT slices based on Slice Location
for ctseri=1:CTSerLen
CTSLoc=zeros(length(CTSeriesFilesInfoCell{ctseri,1}),1);
for filenumber=1:length(CTSeriesFilesInfoCell{ctseri,1})
%     CTSLoc(filenumber,1)=CTSeriesFilesInfoCell{ctseri,1}{filenumber,1}.SliceLocation;
    CTSLoc(filenumber,1)=CTSeriesFilesInfoCell{ctseri,1}{filenumber,1}.ImagePositionPatient(3);
end
[CTSLoc_sorted,CTSLoc_sorted_ind]=sort(CTSLoc,'descend');
CTSeriesFilesInfoCell{ctseri,1}=CTSeriesFilesInfoCell{ctseri,1}(CTSLoc_sorted_ind);
CTSeriesFilesInfoCell{ctseri,1+1}=CTSLoc_sorted;
CTSeriesFilesPathCell{ctseri,1}=CTSeriesFilesPathCell{ctseri,1}(CTSLoc_sorted_ind);
end
end