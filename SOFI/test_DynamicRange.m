clc,clear,close all;
load sofi;
img = sofi{2};
[m,n] = size(img);
kernel = [1,1,1;
          1,1,1;
          1,1,1];
conv_matrix = conv2(img,kernel,'same');
reshape_convMatrix = reshape(conv_matrix,[m*n,1]);                         %�ֳ��о��󣬱��ڷ���
class_matrix = kmeans(reshape_convMatrix,4);                               %k-means�������о���
res_matrix = reshape(class_matrix,[m, n]);                                 %����mask����
figure;
imshow(res_matrix,[]);
disp('Classfication Finished!');
%%        
% end

