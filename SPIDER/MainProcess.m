clc,clear,close all;

image = 'Lab.tif';
[m,n] = size(imread(image));
imgNum = 1000;
imageData = zeros(m,n,imgNum);
for i = 1 : imgNum
    imageData(:,:,i) = imread(image,i);
end

disp('Image infomation was successfully collected!');

[imageStacks, rows, cols] = imgSegment(imageData);
resStacks = cell(m,n);
sparsePos = cell(m,n);
tic;
t1 = cputime;
n = 0;
disp('Sparse Convolution Calculating...');
for i = 1: rows
    for j = 1 : cols
        [resStacks{i,j},sparsePos{i,j}] = spiderMain(imageStacks(i,j), 3);
        n = n + 1;
        fprintf('Progress: %d / %d...... ', n , rows * cols);
        t2 = cputime;
        fprintf('Estimated Finish Time: %s seconds.\n', num2str(((t2-t1)/n) * (rows*cols - n)));
    end
end
fusionImg = cell2mat(resStacks);
imagesp(xy_gray2rgb(fusionImg,'pink'),sprintf('Spider Reconstruction'));

toc;
load chirp;
sound(y,Fs);
