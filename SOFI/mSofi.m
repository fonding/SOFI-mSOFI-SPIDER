function mSofiResult = mSofi(stack,order)
%% MSOFI 
% Moment - SOFI Reconstruction Method.
% 

%% ����������
if nargin == 1 || isempty(order)
    order = 4;                                                             % 4�׾��ع�SOFI = 8���ۻ����ع�SOFI
end

%% �Զ���ȡͼƬ���� ----------����Ŀ�ģ�1.�õ�pic_num����frames��2.�õ�ͼƬ��ά����Ϣ
pic = Tiff(stack);                                                         % TIFF��ջ
frames_info = pic.getTag('ImageDescription');                              % ��ȡͼƬ����Ϣ
newStr = split(frames_info);                                               % ��ȡ��frames����Ϣ
str_splited = split(newStr(2),'=');                                        % �ַ�������
xdim = pic.getTag('ImageLength');                                          % ͼƬ��ά����Ϣ
ydim = pic.getTag('ImageWidth');                                           % ͼƬ��ά����Ϣ
pic_num = str2double(str_splited{2});                                      % TiffͼƬ�а���������ͼ
%% �ɵ��ڲ���
brightness = 3;                                                            % ���ڵ�������ͼƬ�ĶԱȶ�
windowSize = 100;                                                           % local dynamic region compression   %����Window_Size
para.PSF0 = fspecial('gaussian', [4,4], 2);                                % prepare an estimation of the convolution kernel here. In this datasets, cross-correlation doesn't work, therefore we use a rough estimation.
% para.lambda = 2.5;                                                         % DeconvSK parameter
para.lambda = 1.5;                                                         % DeconvSK parameter
para.ItN = 1;                                                              % iteration number.
addpath(genpath(['./Library']))
NF = fspecial('gaussian',[1,21],1);                                        % ����һ���˲���
%% calcualte Moment-information for each block of all frames.
% ImMean_Series = zeros(xdim, ydim);                                       % ����վ���/�洢��Ƭ��ֵ����õ���ͼ��
J = zeros(xdim, ydim, pic_num, 'uint16');                                  % ����ά��Ϊxdim*ydim*200��J��������Ϊuint16
for frameInd = 1 : pic_num
    J(:,:,frameInd) = imread(stack, 'Index', frameInd);                    % ���ζ�ȡÿ��Block��������ͼƬ������
end
% Calculatie average image for this block    
ImMean = mean(double(J), 3);                                               % ��Z���ϼ����ֵ
% ImMean_Series(:, :) = ImMean;                                            % ��ÿ��Block��200֡ͼƬ��ֵ��õ���ͼ����Ϊһ֡���浽ImMean�С�����õ���ͼ�������ɡ�
ImMean_Series(:, :) = ImMean;

    for i = 1: order
%         mSofiResult{i} = momentCalculation(order,brightness,para,windowSize,NF,ImMean,ImMean_Series,xdim,ydim);
        mSofiResult{i} = momentCalculation(order,para,windowSize,NF,ImMean,ImMean_Series,xdim,ydim,pic_num,J);
        disp(['Moment Reconstruction: ',num2str(i),' / ',num2str(order),' Finished!']);
    end
    
end
%% ---------------------------------------------��ֵͼƬ��Ϣ����------------------------------------------------------------

%% ---------------------------------------------����ѭ����ʼ----------------------------------------------------------------
% function momentResult = momentCalculation(order,brightness,para,windowSize,NF,ImMean,ImMean_Series,xdim,ydim)
function momentResult = momentCalculation(order,para,windowSize,NF,ImMean,ImMean_Series,xdim,ydim,pic_num,J)        % brigntnessҪ��Ҫ��ӽ���
%% �߽׾���Ϣ�ļ���
Mseries = zeros(xdim, ydim);                                               % ��ʼ��һ���վ��������洢��Ƭ6�׾��ع�����Ϣ��    
M = zeros(xdim, ydim);                                                     % ����վ��󣬴洢ÿ��Block����6�׾�ʱ�õ�������    
for i0 = 1 : pic_num                                                       % ÿ��Block�к���200֡
    M = M + abs(double(J(:, :, i0)) - ImMean).^sqrt(order);     %�������⣡-�����       % ����ͼ��ĸ߽׾���Ϣ��������������sqrt
%     M = M + (double(J(:, :, i0)) - ImMean).^order;                         %��ͼ��ĸ߽׾���Ϣ��������������sqrt
end
save M;
M = M./pic_num;                                                            % 200��ͼƬ���ֵ
Mseries(:, :) = M;                                                         % ��¼��ÿ��ͼƬ���յľ��ع�����'

%% perform noise filter on M6 along the time axis for each pixel independently.

M_filtered = zeros(xdim, ydim);                                            % ����վ������ڴ洢��Ӧÿ֡ͼƬ�˲����ͼ��
for i0 = 1:xdim                                                            % X�����
    for i1 = 1:ydim                                                        % Y�����
        s = Mseries(i0, i1); s = s(:);                                     % ȡ��ÿ�����ص��Ӧ��20֡ͼƬ����Ϣ��
        s = conv(s,NF,'same');                                             % ������㣬same
        M_filtered(i0, i1) = s;                                            % ���°�s����reshapeΪ1*1*20��ά�ȣ��洢���֮���ͼƬ��
    end
end

%% perform shrinking kernal deconvolution on M6 for each frame.

M_NF_DeconvSK = zeros(xdim, ydim);                                         % ����洢���������Ŀվ���                                                                                                 %��ÿһ֡���б�����
im = M_filtered(:, :);                                                     % ��ȡÿһ֡��ͼƬ��Ϣ��

% now prepare input parameters for xy_DeconvSk

[xdim, ydim] = size(im);                                                   % ��ȡͼƬ��ά��

inputImg = [im, fliplr(im); flipud(im), rot90(im,2)];                      % perform mirror extension to the image in order to surpress ringing artifacts associated with fourier transform due to truncation effect.
para.J0 = inputImg;

output = xy_DeconvSK(para);                                                % ���������
DeconvNFM6 = output(1:xdim,1:ydim);                                        % ѡȡ��ԭͼƬһ����ά�ȡ�
M_NF_DeconvSK(:, :) = DeconvNFM6;                                          % �洢���������֮���ͼƬ��
%% perform then next round of noise filter to the deconvolution result
M_NF_DeconvSK_NF = zeros(xdim, ydim);                                      % ����վ���
for i0 = 1 : xdim                                                          % ����X��
    for i1 = 1 : ydim                                                      % ����Y��
        s = M_NF_DeconvSK(i0,i1,:); s = s(:);                              % ÿ�����ص��ʱ��������Ϣ
        s = conv(s,NF,'same');                                             % ����һ�������˳�����
        M_NF_DeconvSK_NF(i0,i1) = s;                                       % �洢ȥ���������ͼƬ��
    end
end    
%% perform LDRC on the current result and save.
M_NF_DeconvSK_NF_LDRC = zeros(xdim, ydim);                                 % �վ��󣬴洢ʹ��LDRC����֮��Ľ��
                                                                           % ����ÿһ֡ͼ��
InputImage = M_NF_DeconvSK_NF(:,:);                                        % ��ȡÿһ֡ͼƬ����Ϣ                                                                                                       %����LDRC�Ľ���
Mask = ImMean_Series(:,:);                                                 % ��i0��Blockͼ��ľ�ֵ ��ͬ�ϣ�
ldrcResult = xy_QuickLDRC(InputImage, Mask, windowSize);                   % �ó������LDRC��ͼ������
M_NF_DeconvSK_NF_LDRC(:,:) = ldrcResult;                                   % �洢LDRC����֮���ͼƬ
%%
momentResult = M_NF_DeconvSK_NF_LDRC(:, :);                                % �ڶ���ͼƬ��ʾ�߽׾��LDRC����֮���ͼƬ
end                                                                        % ������ֹ

% figure;
% imshow(rot90(mSofiResult*brightness,1));