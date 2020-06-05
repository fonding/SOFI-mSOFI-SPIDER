clc,clear, close all;
load('spider3.mat');
figure;
imshow(fusionImg,[]);
colormap(pink);
result = [];
for i = 1 : 384
    for j = 1 : 384
        if mod(i, 32) == 0 || mod(i, 32) == 1
            result = [result fusionImg(i,j)];
        elseif mod(j, 32) == 0 || mod(j, 32) == 1
            result = [result fusionImg(i,j)];
        end
    end
end
avg = mean(result);
for i = 1 : 384
    for j = 1 : 384
        fusionImg(i,j) = fusionImg(i,j) - avg;
        if fusionImg(i,j) < 0
           fusionImg(i,j) = 0;
        end
    end
end
figure;
imshow(fusionImg,[]);
colormap(pink);