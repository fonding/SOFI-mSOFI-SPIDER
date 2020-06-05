function [imgStack, rows, cols] = imgSegment(img)
[dimension1,dimension2, frames] = size(img);

if dimension1 ~= dimension2
    error('You Must Input A Square Image, (row # == col #) !');
    
elseif dimension1 ~= 128 && dimension1 ~= 64 && dimension1 ~= 256 && dimension1 ~= 1024 && dimension1 ~= 512
    error('Youe input is supposed be a 128*128 or 64*64 or 256*256 or 512 * 512 or 1024*1024 image');
    
else
    count = 0;
    while dimension1 > 32           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dimension1 = dimension1 / 2;
        count = count + 1;
    end
end
if count == 0
    imgStack = cell(1, 1);
    imgStack{1,1} = img;
    rows = dimension1;
    cols = dimension2;
else
    imgStack = cell(2^count,2^count);
%     edgeLength = dimension1 / 2^count;
    rows = 2^count;
    cols = 2^count;
    for i = 1 : rows
        for j = 1 : cols
            imgStack{i,j} = img((i-1)*dimension1+1:i*dimension1,(j-1)*dimension1+1:j*dimension1,:);
        end
    end
end

end