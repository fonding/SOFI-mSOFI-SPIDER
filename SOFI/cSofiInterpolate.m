function csofi = cSofiInterpolate(sofi_image)
%SOFIINTERPOLATE 此处显示有关此函数的摘要
%   此处显示详细说明
% sofi = sofi_image;
[~,n] = size(sofi_image);
[size1,size2] = size(sofi_image{1});
    for i = 1 : n
        csofi{i} = imresize(sofi_image{i},[(size1+3)*i,(size2+3)*i],'bicubic');
    end
end