function res_matrix = classfication(img)
%CLASSFICATION 此处显示有关此函数的摘要
%   此处显示详细说明
% clc,clear,close all;
% load sofi;
% img = sofi{4};
[m,n] = size(img);
kernel = [1,1,1;
          1,1,1;
          1,1,1];
conv_matrix = conv2(img,kernel,'same');
reshape_convMatrix = reshape(conv_matrix,[m*n,1]);                         %分成列矩阵，便于分类
class_matrix = kmeans(reshape_convMatrix,2);                               %k-means方法进行聚类
compare1 = class_matrix == 1;
compare2 = class_matrix == 2;                                              %找到分为第2类的点
if mean(reshape_convMatrix(compare1)) > mean(reshape_convMatrix(compare2)) %对像素点的随机分类进行纠正。使得亮点的分类始终为1，暗点为0
    class_matrix = abs(class_matrix - 2);
else
    class_matrix = abs(class_matrix - 1);
end
res_matrix = reshape(class_matrix,[m, n]);                                 %生成mask矩阵
% figure;
% imshow(res_matrix);
disp('Classfication Finished!');
%%        
% se = strel('sphere',13);
% closeBW = imclose(res_matrix,se);
% figure, imshow(closeBW)
end

