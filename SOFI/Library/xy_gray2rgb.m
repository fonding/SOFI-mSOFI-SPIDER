% This code is distributed under MIT license, please refer to the LICENSE file in the package for details.
%
%   Copyright (c) 2018 Xiyu Yi
%
%   Author of the code: Xiyu Yi
%   Email of the author: xiyu.yi@gmail.com

function output = xy_gray2rgb(im, map)
    if nargin < 2 || isempty(map)
        map = 'gray';
    end
    if ischar(map)
        map = [0,0,0;imresize(colormap(eval(map)), [256, 3])];                         % 伪彩，在colormap第一行加一行(0,0,0)
    else
        error('The 2nd parameter of xy_gray2rgb function must be a string!');
    end
im(isnan(im))= eps;
im = double(im);
im = (im - min(im(:)))./(max(im(:))-min(im(:))).*length(map(:,1));
im = round(im);
output = ind2rgb(im, map);
