function res_matrix = classfication(img)
%CLASSFICATION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% clc,clear,close all;
% load sofi;
% img = sofi{4};
[m,n] = size(img);
kernel = [1,1,1;
          1,1,1;
          1,1,1];
conv_matrix = conv2(img,kernel,'same');
reshape_convMatrix = reshape(conv_matrix,[m*n,1]);                         %�ֳ��о��󣬱��ڷ���
class_matrix = kmeans(reshape_convMatrix,2);                               %k-means�������о���
compare1 = class_matrix == 1;
compare2 = class_matrix == 2;                                              %�ҵ���Ϊ��2��ĵ�
if mean(reshape_convMatrix(compare1)) > mean(reshape_convMatrix(compare2)) %�����ص�����������о�����ʹ������ķ���ʼ��Ϊ1������Ϊ0
    class_matrix = abs(class_matrix - 2);
else
    class_matrix = abs(class_matrix - 1);
end
res_matrix = reshape(class_matrix,[m, n]);                                 %����mask����
% figure;
% imshow(res_matrix);
disp('Classfication Finished!');
%%        
% se = strel('sphere',13);
% closeBW = imclose(res_matrix,se);
% figure, imshow(closeBW)
end

