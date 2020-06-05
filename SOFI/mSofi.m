function mSofiResult = mSofi(stack,order)
%% MSOFI 
% Moment - SOFI Reconstruction Method.
% 

%% 检查输入参数
if nargin == 1 || isempty(order)
    order = 4;                                                             % 4阶距重构SOFI = 8阶累积量重构SOFI
end

%% 自动获取图片参数 ----------最终目的：1.得到pic_num，即frames。2.得到图片的维度信息
pic = Tiff(stack);                                                         % TIFF堆栈
frames_info = pic.getTag('ImageDescription');                              % 获取图片的信息
newStr = split(frames_info);                                               % 提取出frames的信息
str_splited = split(newStr(2),'=');                                        % 字符串分离
xdim = pic.getTag('ImageLength');                                          % 图片的维度信息
ydim = pic.getTag('ImageWidth');                                           % 图片的维度信息
pic_num = str2double(str_splited{2});                                      % Tiff图片中包含多少张图
%% 可调节参数
brightness = 3;                                                            % 用于调节生成图片的对比度
windowSize = 100;                                                           % local dynamic region compression   %定义Window_Size
para.PSF0 = fspecial('gaussian', [4,4], 2);                                % prepare an estimation of the convolution kernel here. In this datasets, cross-correlation doesn't work, therefore we use a rough estimation.
% para.lambda = 2.5;                                                         % DeconvSK parameter
para.lambda = 1.5;                                                         % DeconvSK parameter
para.ItN = 1;                                                              % iteration number.
addpath(genpath(['./Library']))
NF = fspecial('gaussian',[1,21],1);                                        % 定义一个滤波器
%% calcualte Moment-information for each block of all frames.
% ImMean_Series = zeros(xdim, ydim);                                       % 定义空矩阵/存储照片均值处理得到的图像。
J = zeros(xdim, ydim, pic_num, 'uint16');                                  % 定义维度为xdim*ydim*200的J矩阵，类型为uint16
for frameInd = 1 : pic_num
    J(:,:,frameInd) = imread(stack, 'Index', frameInd);                    % 依次读取每个Block里面所有图片的数据
end
% Calculatie average image for this block    
ImMean = mean(double(J), 3);                                               % 在Z轴上计算均值
% ImMean_Series(:, :) = ImMean;                                            % 把每个Block中200帧图片均值后得到的图像作为一帧，存到ImMean中。处理得到的图像计算完成。
ImMean_Series(:, :) = ImMean;

    for i = 1: order
%         mSofiResult{i} = momentCalculation(order,brightness,para,windowSize,NF,ImMean,ImMean_Series,xdim,ydim);
        mSofiResult{i} = momentCalculation(order,para,windowSize,NF,ImMean,ImMean_Series,xdim,ydim,pic_num,J);
        disp(['Moment Reconstruction: ',num2str(i),' / ',num2str(order),' Finished!']);
    end
    
end
%% ---------------------------------------------均值图片信息结束------------------------------------------------------------

%% ---------------------------------------------函数循环开始----------------------------------------------------------------
% function momentResult = momentCalculation(order,brightness,para,windowSize,NF,ImMean,ImMean_Series,xdim,ydim)
function momentResult = momentCalculation(order,para,windowSize,NF,ImMean,ImMean_Series,xdim,ydim,pic_num,J)        % brigntness要不要添加进来
%% 高阶距信息的计算
Mseries = zeros(xdim, ydim);                                               % 初始化一个空矩阵，用来存储照片6阶距重构的信息。    
M = zeros(xdim, ydim);                                                     % 定义空矩阵，存储每个Block计算6阶距时得到的数据    
for i0 = 1 : pic_num                                                       % 每个Block中含有200帧
    M = M + abs(double(J(:, :, i0)) - ImMean).^sqrt(order);     %复数问题！-待解决       % 计算图像的高阶距信息。！！！！！！sqrt
%     M = M + (double(J(:, :, i0)) - ImMean).^order;                         %算图像的高阶距信息。！！！！！！sqrt
end
save M;
M = M./pic_num;                                                            % 200张图片求均值
Mseries(:, :) = M;                                                         % 记录下每张图片最终的距重构数据'

%% perform noise filter on M6 along the time axis for each pixel independently.

M_filtered = zeros(xdim, ydim);                                            % 定义空矩阵，用于存储对应每帧图片滤波后的图像
for i0 = 1:xdim                                                            % X轴遍历
    for i1 = 1:ydim                                                        % Y轴遍历
        s = Mseries(i0, i1); s = s(:);                                     % 取出每个像素点对应的20帧图片的信息。
        s = conv(s,NF,'same');                                             % 卷积计算，same
        M_filtered(i0, i1) = s;                                            % 重新把s矩阵reshape为1*1*20的维度，存储卷积之后的图片。
    end
end

%% perform shrinking kernal deconvolution on M6 for each frame.

M_NF_DeconvSK = zeros(xdim, ydim);                                         % 定义存储反卷积计算的空矩阵。                                                                                                 %对每一帧进行遍历。
im = M_filtered(:, :);                                                     % 读取每一帧的图片信息。

% now prepare input parameters for xy_DeconvSk

[xdim, ydim] = size(im);                                                   % 获取图片的维度

inputImg = [im, fliplr(im); flipud(im), rot90(im,2)];                      % perform mirror extension to the image in order to surpress ringing artifacts associated with fourier transform due to truncation effect.
para.J0 = inputImg;

output = xy_DeconvSK(para);                                                % 反卷积计算
DeconvNFM6 = output(1:xdim,1:ydim);                                        % 选取与原图片一样的维度。
M_NF_DeconvSK(:, :) = DeconvNFM6;                                          % 存储反卷积计算之后的图片。
%% perform then next round of noise filter to the deconvolution result
M_NF_DeconvSK_NF = zeros(xdim, ydim);                                      % 定义空矩阵
for i0 = 1 : xdim                                                          % 遍历X轴
    for i1 = 1 : ydim                                                      % 遍历Y轴
        s = M_NF_DeconvSK(i0,i1,:); s = s(:);                              % 每个像素点的时间序列信息
        s = conv(s,NF,'same');                                             % 再做一遍卷积，滤除噪声
        M_NF_DeconvSK_NF(i0,i1) = s;                                       % 存储去除噪声后的图片。
    end
end    
%% perform LDRC on the current result and save.
M_NF_DeconvSK_NF_LDRC = zeros(xdim, ydim);                                 % 空矩阵，存储使用LDRC计算之后的结果
                                                                           % 遍历每一帧图像
InputImage = M_NF_DeconvSK_NF(:,:);                                        % 提取每一帧图片的信息                                                                                                       %定义LDRC的阶数
Mask = ImMean_Series(:,:);                                                 % 第i0个Block图像的均值 （同上）
ldrcResult = xy_QuickLDRC(InputImage, Mask, windowSize);                   % 得出计算的LDRC的图像结果。
M_NF_DeconvSK_NF_LDRC(:,:) = ldrcResult;                                   % 存储LDRC计算之后的图片
%%
momentResult = M_NF_DeconvSK_NF_LDRC(:, :);                                % 第二张图片显示高阶距和LDRC计算之后的图片
end                                                                        % 函数终止

% figure;
% imshow(rot90(mSofiResult*brightness,1));