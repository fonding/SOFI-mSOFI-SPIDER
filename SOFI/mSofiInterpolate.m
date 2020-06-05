function msofi = mSofiInterpolate(sofi_image)
%SOFIINTERPOLATE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[~,n] = size(sofi_image);
[size1,size2] = size(sofi_image{1});
    for i = 1 : n
        msofi{i} = imresize(sofi_image{i},[(size1)*i,(size2)*i],'bicubic');
    end
end