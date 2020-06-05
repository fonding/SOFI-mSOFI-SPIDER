function csofi = cSofiInterpolate(sofi_image)
%SOFIINTERPOLATE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% sofi = sofi_image;
[~,n] = size(sofi_image);
[size1,size2] = size(sofi_image{1});
    for i = 1 : n
        csofi{i} = imresize(sofi_image{i},[(size1+3)*i,(size2+3)*i],'bicubic');
    end
end