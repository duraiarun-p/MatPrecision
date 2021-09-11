clc;clear;close all;
load('/home/arun/Documents/MATLAB/ImageDB/PrintoutDB/ZC001-DB2.mat','PlanCTData','PosCBCTData');
%%
CT1=PlanCTData(:,:,200);
CT2=PosCBCTData(:,:,44);
% CT=imabsdiff(CT1,CT2);
% CT1l=log(CT1);
% CT2l=log(CT2);
% CTl=imabsdiff(CT1l,CT2l);
% %%
% figure(1);subplot(131),imshow(CT1,[]);subplot(132),imshow(CT2,[]);subplot(133),imshow(CT,[]);
% figure(2);subplot(131),imshow(CT1l,[]);subplot(132),imshow(CT2l,[]);subplot(133),imshow(CTl,[]);
%%
inputlayer = imageInputLayer([512 512 1],'Name','input');
layer1 = convolution2dLayer(3,128,'Padding','same');
relulayer = reluLayer('Name','relu1');
maxpoollayer = averagePooling2dLayer(1,'Padding','same');

layers = [ 
          inputlayer
          layer1
          relulayer
          maxpoollayer
          fullyConnectedLayer(10)
          softmaxLayer
          classificationLayer];
      
% analyzeNetwork(layers);
options = trainingOptions('sgdm', ...
    'MaxEpochs',5,...
    'InitialLearnRate',1e-4, ...
    'Verbose',false, ...
    'Plots','training-progress');
%%
XTrain=cat(3,PlanCTData(:,:,51:60),PosCBCTData(:,:,44));
YTrain=[ones(5,1);zeros(5,1)];

%%
% net = trainNetwork(XTrain,YTrain,layers,options);