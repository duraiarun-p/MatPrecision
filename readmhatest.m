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

filename='/home/arun/Documents/MATLAB/ImageOutputs/PlastiOutputs/demons_vf.mha';

%%
% function data = mhd_read_image(filename)
% info = mha_read_header(filename);
% data = mha_read_volume(info);
% end