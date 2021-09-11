clc;clear;close all;
%%
cd /home/arun/Documents/PyWSPrecision/dblm
% [MTE, metaM] = nrrdread('MTEimg.nrrd');
% [ITE, metaI] = nrrdread('ITEimg.nrrd');
load pydblmout.mat
Dice=mydice(MTE,MT2);


%%
function DiceMatlab=mydice(NMT,IOT)
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
end