clc,clear,close all;

addpath(genpath(['./Library']))
NF = fspecial('gaussian',[1,21],2);                                                     %����һ���˲���
% get dimensions
im = imread(['./Demo_dataset/Block1.tif']);                                             %��ȡͼƬ
[xdim, ydim] = size(im);                                                                %��ȡͼƬ��dimension��Ϣ

%% calcualte M6 for each block of 200 frames, and 20 blocks in total.
M6series = zeros(xdim, ydim, 20);                                                       %��ʼ��һ���վ��������洢��Ƭ6�׾��ع�����Ϣ��
ImMean_Series = zeros(xdim, ydim, 20);                                                  %����վ���/�洢��Ƭ��ֵ����õ���ͼ��

for blockN = 1 : 20
    J = zeros(xdim, ydim, 200, 'uint16');                                               %����ά��Ϊxdim*ydim*200��J��������Ϊuint16
    disp(['calculating M6 of block ',num2str(blockN),'/20...'])                         %��ʾ������
    for frameInd = 1 : 200
        im = imread(['./Demo_dataset/Block',num2str(blockN),'.tif'], 'Index', frameInd);%���ζ�ȡÿ��Block�����200��ͼƬ����
        J(:,:,frameInd) = im;                                                           %�����ݴ浽J������
    end                                                                              
    %   calculatie average image for this block
    ImMean = mean(double(J), 3);                                                        %��Z���ϼ����ֵ
    ImMean_Series(:, :, blockN) = ImMean;                                               %��ÿ��Block��200֡ͼƬ��ֵ��õ���ͼ����Ϊһ֡���浽ImMean�С�����õ���ͼ�������ɡ�
    
    %   calculate M6 of this block
    M6 = zeros(xdim, ydim);                                                             %����վ��󣬴洢ÿ��Block����6�׾�ʱ�õ�������
    
    for i0 = 1 : 200                                                                    %ÿ��Block�к���200֡
        M6 = M6 + (double(J(:, :, i0)) - ImMean).^6; % directly compute.               %����ͼ���6�׾�
    end
    M6 = M6./200;                                                                       %200��ͼƬ���ֵ
    M6series(:, :, blockN) = M6;                                                        %��¼��ÿ��ͼƬ���յľ��ع�����
end

%% perform noise filter on M6 along the time axis for each pixel independently.
M6_filtered = zeros(xdim, ydim, 20);                                                                            %����վ������ڴ洢��Ӧÿ֡ͼƬ�˲����ͼ��
for i0 = 1:xdim                                                                                                 %X�����
    if rem(i0, 100) == 0                                                                                        %rem(a,b),������b����a֮�������
    disp(['Apply noise filter on M6 series, row (',num2str(i0),'-',num2str(min(i0+100-1,283)),')/283...'])      %��ʾ������
    end                                                                                                         %
    for i1 = 1:ydim                                                                                             %Y�����
        s = M6series(i0, i1, :); s = s(:); % take the time sequence                                             %ȡ��ÿ�����ص��Ӧ��20֡ͼƬ����Ϣ��
        s = conv(s,NF,'same'); % filter.                                                                        %������㣬same
        M6_filtered(i0, i1, :) = reshape(s, [1, 1, 20]); % store into the filtered variable                     %���°�s����reshapeΪ1*1*20��ά�ȣ��洢���֮���ͼƬ��
    end
end

%% perform shrinking kernal deconvolution on M6 for each frame.
M6_NF_DeconvSK = zeros(xdim, ydim, 20);                                                                         %����洢���������Ŀվ���
for i0 = 1 : 20                                                                                                 %��ÿһ֡���б�����
    im = M6_filtered(:, :, i0);                                                                                 %��ȡÿһ֡��ͼƬ��Ϣ��
    % now prepare input parameters for xy_DeconvSk
        [xdim, ydim] = size(im);                                                                                %��ȡͼƬ��ά��
        inputImg = [im, fliplr(im); flipud(im), rot90(im,2)];% perform mirror extension to the image in order to surpress ringing artifacts associated with fourier transform due to truncation effect.
        para.J0 = inputImg;
        para.PSF0 = fspecial('gaussian', [51,51], 2); % prepare an estimation of the convolution kernel here. In this datasets, cross-correlation doesn't work, therefore we use a rough estimation.
        para.lambda = 1.5; %DeconvSK parameter
        para.ItN = 20; % iteration number.
    disp(['calculating DeconvSK on block ',num2str(i0),'/20...'])                                               %��ʾ������
    output = xy_DeconvSK(para);                                                                                 %���������
    DeconvNFM6 = output(1:xdim, 1:ydim);                                                                        %ѡȡ��ԭͼƬһ����ά�ȡ�
    M6_NF_DeconvSK(:, :, i0) = DeconvNFM6;                                                                      %�洢���������֮���ͼƬ��
end
%% perform then next round of noise filter to the deconvolution result
M6_NF_DeconvSK_NF = zeros(xdim, ydim, 20);                                                                            %����վ���
for i0 = 1 : xdim                                                                                                     %����X��
    if rem(i0, 100) == 0                                                                                              %��ʾ�������
    disp(['Apply noise filter on deconvolved M6 series, row (',num2str(i0),'-',num2str(min(i0+100-1,283)),')/283...'])%
    end
    for i1 = 1 : ydim                                                                                                 %����Y��
        s = M6_NF_DeconvSK(i0,i1,:); s = s(:);                                                                        %ÿ�����ص��ʱ��������Ϣ
        s = conv(s,NF,'same');                                                                                        %����һ�������˳�����
        M6_NF_DeconvSK_NF(i0,i1,:) = reshape(s,[1,1,20]);                                                             %�洢ȥ���������ͼƬ��
    end
end    
%% perform LDRC on the current result and save.
M6_NF_DeconvSK_NF_LDRC = zeros(xdim, ydim, 20);                                                                       %�վ��󣬴洢ʹ��LDRC����֮��Ľ��
for i0 = 1 : 20                                                                                                       %����ÿһ֡ͼ��
    InputImage = M6_NF_DeconvSK_NF(:,:,i0);                                                                           %��ȡÿһ֡ͼƬ����Ϣ
    order = 7;                                                                                                        %����LDRC�Ľ���
    Mask = ImMean_Series(:,:,i0);                                                                                     %��i0��Blockͼ��ľ�ֵ ��ͬ�ϣ�
    windowSize = 25;                                                                                                  %����Window_Size
    disp(['calculate LDRC for frame ',num2str(i0),'/20'])                                                             %��ʾ���������
    ldrcResult = xy_QuickLDRC(InputImage, Mask, windowSize);                                                          %�ó������LDRC��ͼ������
    M6_NF_DeconvSK_NF_LDRC(:,:,i0) = ldrcResult;                                                                      %�洢LDRC����֮���ͼƬ
end

%% save the result
save demo_result.mat M6_NF_DeconvSK_NF_LDRC ImMean_Series                                                             %�洢���

%% produce visualization of the result        20֡ͼ������һ��video���ӻ���
vid = VideoWriter('12_order_mSOFI.avi');
vid.Quality = 100;                                                                                                      %��Ƶ��������
vid.FrameRate=10;                                                                                                       %��Ƶ���ŵ�֡��

open(vid);                                                                                                              %�½�һ��video
cmap = colormap(pink);                                                                                                  %α��
cmap = [0,0,0;imresize(cmap, [256, 3])];                                                                                %��colormap��һ�м�һ��(0,0,0)
for i0 = 1 : 20
    im1 = xy_gray2rgb(ImMean_Series(:, :, i0), cmap);                                                                   %��һ��ͼƬ��ʾ��ֵͼƬ��
    im2 = xy_gray2rgb(M6_NF_DeconvSK_NF_LDRC(:, :, i0), cmap);                                                          %�ڶ���ͼƬ��ʾ�߽׾��LDRC����֮���ͼƬ
    
    figure(1); imshow([im1, im2.* 1.4]);%adjust display contrast for better visualization                                %ͬʱ��ʾ������ͼƬ����ֵͼƬ����ߣ����ع�ͼƬ���ұ�
    drawnow;                                                                                                            %����ͼ��
    c = getimage;                                                                                                       %����ͼ��
    writeVideo(vid, c);                                                                                                 %��ͼ��cд���Vid��Video����
end
close(vid);                                                                                                              %END
