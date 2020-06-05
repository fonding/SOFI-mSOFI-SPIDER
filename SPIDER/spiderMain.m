function [imgSegments, SparsePositions] = spiderMain(img, zoom)
%SPIDERMAIN 此处显示有关此函数的摘要
%   此处显示详细说明
% *****************************************************************************************************************
% SParse deconvolution of hIgh-Density supER-resolution images
% 建议所有需要修改的参数都在这个文件里修改
% *****************************************************************************************************************
%% Load the data
data = cell2mat(img);
[dimension1, dimension2, frames] = size(data);
frames = 100;
% dataf = zeros(dimension1,dimension2,frames);
kappa = zeros(frames,1);
for i = 1: frames
%     kappa(i) = 1/mean(mean(data(:,:,i)));
    kappa(i) = mean(mean(data(:,:,i)));
end
%% Description
% Initialisation of the parameters needed for every map.
% As this demoset doesn't require the images to be broken up into different patches, 
% there is no need for an overlap to compensate for emitters at
% the border of every patch. This also doesn't require to have a subvector 
% that says which patches should be calculated and which ones can be 
% skipped. The data has also been simulated so no baseline removal is needed.
%% Perform SPIDER calculations
sigma = 1;                                                                 % Sigma has been chosen so that it's equal to the size of a pixel
% zoom = 4;                                                                  % Oversampling factor (usually between 1 - 8)
overlap = 0;                                                               % No need of an overlap as we have no patches
% PatchSize = 30;                                                            % Size of every patch (here it's the size of the image)
PatchSize = dimension1;

BaselineRemoval = 0;                                                       % No need to remove a baseline
% kappa = [1100; 950; 850; 700; 650; 550; 550; 500; 450; 400; 300; 250; 300]; % Different densities require different penalty parameters
subvector = 1;                                                             % No patches require a vector with 1
SparseImage = zeros(dimension1*zoom, dimension2*zoom, frames);
SparsePositions = cell(1,frames);
% parfor i = 1:frames
parfor i = 1 : 10               %这个的大小影响了图像的连续性。i为多少代表用多少张图象重建
    [SparseImage(:,:,i), SparsePositions{i}] = SpiderInterface(data(:,:,i), sigma, zoom, overlap, PatchSize, BaselineRemoval, kappa(i,:), subvector);
end
imgSegments = sum(SparseImage,3);
end