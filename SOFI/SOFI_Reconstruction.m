%% Image TIFF file or image sequence.

clc,clear,close all;
stack='Lab.tif';

%% Analyze the image sequence.

[sofi,grid]=sofiCumulants(stack);
[sofi,fwhm]=sofiFlatten([],sofi,grid);
[ratio,density,brightness]=sofiParameters(sofi);   % Need flat cumulants to get parameters
sofi=sofiLinearize(sofi,fwhm);                     % "Blind" linearization (no parameters)
sofi=cSofiInterpolate(sofi);                       % "Interpolate" the image to return the original size
m_Sofi_Result = mSofi(stack,8);                    % Stack is the image's name, the 2nd input is the moment-order;
mSofi = mSofiInterpolate(m_Sofi_Result);           % "interpolate" the image to return the interpolated size
for i = 1 : 8
    if i == 1
        img{1} = sofi{1};
    else
        img{i}=sofiBalance(sofi,ratio,sofi{i-1},sofi{i});       % b-sofi¶ÑÕ»
    end
end
%% Display the results.
for n=1:8
   imagesp(xy_gray2rgb(img{n},'pink'),sprintf('bSOFI - %d. order',n));
   imagesp(xy_gray2rgb(sofi{n}.*classfication(img{n}),'pink'),sprintf('Modified Cumulant SOFI - %d. order',n));
   imagesp(xy_gray2rgb(sofi{n},'pink'),sprintf('Cumulant SOFI - %d. order',n));
   imagesp(xy_gray2rgb(mSofi{n},'pink'),sprintf('Moment SOFI - %d. order',n));
   imagesp(xy_gray2rgb(mSofi{n}.*classfication(img{n}),'pink'),sprintf('Modified Moment SOFI - %d. order',n));
end
% imagesp(img,'Balanced');
% figure;
disp('Finished!');
load chirp;
sound(y,Fs);

