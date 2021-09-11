clc;clear;close all;
% profile on;
addpath('/home/arun/Documents/MATLAB/MAWS_Precision');
addpath('/home/arun/Documents/MATLAB/DicomGUI');
addpath('/home/arun/Documents/MATLAB/DicomGUI/Scripting');

mypath='/home/arun/Documents/MATLAB/ImageDB/PrintoutDB/';
cd(mypath);
%%
tic;
load('testprintout008.mat');
timeL=toc;