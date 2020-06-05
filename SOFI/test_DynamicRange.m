clc,clear,close all;
load sofi;
img = sofi{2};
[m,n] = size(img);
kernel = [1,1,1;
          1,1,1;
          1,1,1];
conv_matrix = conv2(img,kernel,'same');
reshape_convMatrix = reshape(conv_matrix,[m*n,1]);                         %分成列矩阵，便于分类
class_matrix = kmeans(reshape_convMatrix,4);                               %k-means方法进行聚类
res_matrix = reshape(class_matrix,[m, n]);                                 %生成mask矩阵
figure;
imshow(res_matrix,[]);
disp('Classfication Finished!');
%%        
% end

