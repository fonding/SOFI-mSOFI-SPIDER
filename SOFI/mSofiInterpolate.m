function msofi = mSofiInterpolate(sofi_image)
%SOFIINTERPOLATE 此处显示有关此函数的摘要
%   此处显示详细说明
[~,n] = size(sofi_image);
[size1,size2] = size(sofi_image{1});
    for i = 1 : n
        msofi{i} = imresize(sofi_image{i},[(size1)*i,(size2)*i],'bicubic');
    end
end