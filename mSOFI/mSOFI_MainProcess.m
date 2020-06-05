clc,clear,close all;

addpath(genpath(['./Library']))
NF = fspecial('gaussian',[1,21],2);                                                     %定义一个滤波器
% get dimensions
im = imread(['./Demo_dataset/Block1.tif']);                                             %读取图片
[xdim, ydim] = size(im);                                                                %获取图片的dimension信息

%% calcualte M6 for each block of 200 frames, and 20 blocks in total.
M6series = zeros(xdim, ydim, 20);                                                       %初始化一个空矩阵，用来存储照片6阶距重构的信息。
ImMean_Series = zeros(xdim, ydim, 20);                                                  %定义空矩阵/存储照片均值处理得到的图像。

for blockN = 1 : 20
    J = zeros(xdim, ydim, 200, 'uint16');                                               %定义维度为xdim*ydim*200的J矩阵，类型为uint16
    disp(['calculating M6 of block ',num2str(blockN),'/20...'])                         %显示进度条
    for frameInd = 1 : 200
        im = imread(['./Demo_dataset/Block',num2str(blockN),'.tif'], 'Index', frameInd);%依次读取每个Block里面的200张图片数据
        J(:,:,frameInd) = im;                                                           %将数据存到J矩阵中
    end                                                                              
    %   calculatie average image for this block
    ImMean = mean(double(J), 3);                                                        %在Z轴上计算均值
    ImMean_Series(:, :, blockN) = ImMean;                                               %把每个Block中200帧图片均值后得到的图像作为一帧，存到ImMean中。处理得到的图像计算完成。
    
    %   calculate M6 of this block
    M6 = zeros(xdim, ydim);                                                             %定义空矩阵，存储每个Block计算6阶距时得到的数据
    
    for i0 = 1 : 200                                                                    %每个Block中含有200帧
        M6 = M6 + (double(J(:, :, i0)) - ImMean).^6; % directly compute.               %计算图像的6阶距
    end
    M6 = M6./200;                                                                       %200张图片求均值
    M6series(:, :, blockN) = M6;                                                        %记录下每张图片最终的距重构数据
end

%% perform noise filter on M6 along the time axis for each pixel independently.
M6_filtered = zeros(xdim, ydim, 20);                                                                            %定义空矩阵，用于存储对应每帧图片滤波后的图像
for i0 = 1:xdim                                                                                                 %X轴遍历
    if rem(i0, 100) == 0                                                                                        %rem(a,b),返回用b除以a之后的余数
    disp(['Apply noise filter on M6 series, row (',num2str(i0),'-',num2str(min(i0+100-1,283)),')/283...'])      %显示进度条
    end                                                                                                         %
    for i1 = 1:ydim                                                                                             %Y轴遍历
        s = M6series(i0, i1, :); s = s(:); % take the time sequence                                             %取出每个像素点对应的20帧图片的信息。
        s = conv(s,NF,'same'); % filter.                                                                        %卷积计算，same
        M6_filtered(i0, i1, :) = reshape(s, [1, 1, 20]); % store into the filtered variable                     %重新把s矩阵reshape为1*1*20的维度，存储卷积之后的图片。
    end
end

%% perform shrinking kernal deconvolution on M6 for each frame.
M6_NF_DeconvSK = zeros(xdim, ydim, 20);                                                                         %定义存储反卷积计算的空矩阵。
for i0 = 1 : 20                                                                                                 %对每一帧进行遍历。
    im = M6_filtered(:, :, i0);                                                                                 %读取每一帧的图片信息。
    % now prepare input parameters for xy_DeconvSk
        [xdim, ydim] = size(im);                                                                                %获取图片的维度
        inputImg = [im, fliplr(im); flipud(im), rot90(im,2)];% perform mirror extension to the image in order to surpress ringing artifacts associated with fourier transform due to truncation effect.
        para.J0 = inputImg;
        para.PSF0 = fspecial('gaussian', [51,51], 2); % prepare an estimation of the convolution kernel here. In this datasets, cross-correlation doesn't work, therefore we use a rough estimation.
        para.lambda = 1.5; %DeconvSK parameter
        para.ItN = 20; % iteration number.
    disp(['calculating DeconvSK on block ',num2str(i0),'/20...'])                                               %显示进度条
    output = xy_DeconvSK(para);                                                                                 %反卷积计算
    DeconvNFM6 = output(1:xdim, 1:ydim);                                                                        %选取与原图片一样的维度。
    M6_NF_DeconvSK(:, :, i0) = DeconvNFM6;                                                                      %存储反卷积计算之后的图片。
end
%% perform then next round of noise filter to the deconvolution result
M6_NF_DeconvSK_NF = zeros(xdim, ydim, 20);                                                                            %定义空矩阵
for i0 = 1 : xdim                                                                                                     %遍历X轴
    if rem(i0, 100) == 0                                                                                              %显示计算进度
    disp(['Apply noise filter on deconvolved M6 series, row (',num2str(i0),'-',num2str(min(i0+100-1,283)),')/283...'])%
    end
    for i1 = 1 : ydim                                                                                                 %遍历Y轴
        s = M6_NF_DeconvSK(i0,i1,:); s = s(:);                                                                        %每个像素点的时间序列信息
        s = conv(s,NF,'same');                                                                                        %再做一遍卷积，滤除噪声
        M6_NF_DeconvSK_NF(i0,i1,:) = reshape(s,[1,1,20]);                                                             %存储去除噪声后的图片。
    end
end    
%% perform LDRC on the current result and save.
M6_NF_DeconvSK_NF_LDRC = zeros(xdim, ydim, 20);                                                                       %空矩阵，存储使用LDRC计算之后的结果
for i0 = 1 : 20                                                                                                       %遍历每一帧图像
    InputImage = M6_NF_DeconvSK_NF(:,:,i0);                                                                           %提取每一帧图片的信息
    order = 7;                                                                                                        %定义LDRC的阶数
    Mask = ImMean_Series(:,:,i0);                                                                                     %第i0个Block图像的均值 （同上）
    windowSize = 25;                                                                                                  %定义Window_Size
    disp(['calculate LDRC for frame ',num2str(i0),'/20'])                                                             %显示计算进度条
    ldrcResult = xy_QuickLDRC(InputImage, Mask, windowSize);                                                          %得出计算的LDRC的图像结果。
    M6_NF_DeconvSK_NF_LDRC(:,:,i0) = ldrcResult;                                                                      %存储LDRC计算之后的图片
end

%% save the result
save demo_result.mat M6_NF_DeconvSK_NF_LDRC ImMean_Series                                                             %存储结果

%% produce visualization of the result        20帧图像，做成一个video可视化。
vid = VideoWriter('12_order_mSOFI.avi');
vid.Quality = 100;                                                                                                      %视频质量参数
vid.FrameRate=10;                                                                                                       %视频播放的帧率

open(vid);                                                                                                              %新建一个video
cmap = colormap(pink);                                                                                                  %伪彩
cmap = [0,0,0;imresize(cmap, [256, 3])];                                                                                %在colormap第一行加一行(0,0,0)
for i0 = 1 : 20
    im1 = xy_gray2rgb(ImMean_Series(:, :, i0), cmap);                                                                   %第一张图片显示均值图片。
    im2 = xy_gray2rgb(M6_NF_DeconvSK_NF_LDRC(:, :, i0), cmap);                                                          %第二张图片显示高阶距和LDRC计算之后的图片
    
    figure(1); imshow([im1, im2.* 1.4]);%adjust display contrast for better visualization                                %同时显示出两张图片，均值图片在左边，距重构图片在右边
    drawnow;                                                                                                            %更新图窗
    c = getimage;                                                                                                       %捕获图像
    writeVideo(vid, c);                                                                                                 %将图窗c写入进Vid（Video）里
end
close(vid);                                                                                                              %END
