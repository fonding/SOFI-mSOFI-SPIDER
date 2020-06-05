% clc,clear,close all;
% load('sofi.mat');
% load('spider3.mat');
clc,clear,close all;
load('spider4.mat');
load('sofi.mat');
sofi8 = img{8};
spider4 = fusionImg;
sofi4 = img{4};

result = [];
for i = 1 : 512
    for j = 1 : 512
        if mod(i, 32) == 0 || mod(i, 32) == 1
            result = [result spider4(i,j)];
        elseif mod(j, 32) == 0 || mod(j, 32) == 1
            result = [result spider4(i,j)];
        end
    end
end
avg = mean(result);
for i = 1 : 512
    for j = 1 : 512
        spider4(i,j) = spider4(i,j) - avg;
        if spider4(i,j) < 0
           spider4(i,j) = 0;
        end
    end
end
% figure;
% imshow(spider4,[]);
% colormap(pink);

figure;
imshow(sofi4,[]);
%% x = 256 y = 660      ------------ x = 400 y = 660
spider4 = normalize(spider4,'range',[0, 1]);
sofi8 = normalize(sofi8,'range',[0, 1]);
sofi4 = normalize(sofi4,'range',[0, 1]);
sofi8 = imresize(sofi8,[1023,1023]);
sofi4 = interp2(sofi4,1);
spider4 = interp2(spider4,1);

spider_4 = spider4(256:400, 661:661);
sofi_4 = sofi4(256:400, 661:661);
sofi_8 = sofi8(256:400, 661:661);
figure;
subplot(3,1,1);
bar(spider_4);
subplot(3,1,2);
bar(sofi_4);
subplot(3,1,3);
bar(sofi_8);
