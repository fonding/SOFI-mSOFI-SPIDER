% *************************************************************************
% SParse deconvolution of hIgh-Density supER-resolution images 
% *************************************************************************
%
% Reference to the publication:
%   Hugelier, S., de Rooi, J.J., Bernex, R., Duwé, S., Devos, O., Sliwa, M.,
%   Dedecker, P., Eilers, P.H.C. & Ruckebusch, C. (2016). Sparse 
%   deconvolution of high-density super-resolution images, Sci. Rep. 6, 
%   21413; doi: 10.1038/srep21413 (2016)
%
% Program to calculate four figures of merit of a fit (for simulated
% data), which are Recall (in %), Accuracy (in nm), the number of False
% Positives (in arbitrary units) and the Sparsity (in %). 
%
% Authors:  
%       - HUGELIER, S.      (1)
%       - DEVOS, O.         (1)
%       - RUCKEBUSCH, C.    (1)
%           
% (1):  LAboratoire de Spectrochimie Infrarouge et Raman (LASIR)
%       Université de Lille 1, UMR CNRS 8516
%       Bât C5, Cité Scientifique
%       59655 Villeneuve d'Ascq - France
% *************************************************************************

function [Recall, Accuracy, FalsePositive, Sparsity, MeanCentered] = FOM_Test(position, spider, PixelToNm, Acc, Zoom)

% *************************************************************************
% *************************************************************************
% 
% Program to test the Spider method. It calculated the mean centered 
% emitters first and it then allows the user to determine some 
% characteristics to determine whether the program is calculating the 
% correct things or not.
% 
% Input:
% 
%   Position:       The matrix with the simulated positions of the emitters
%   Spider:         The results coming from Spider (Positions table in the
%                   Spider routine)
%   PixelToNm:      The conversion of pixels to nm (e.g.: 1 pixel = 150nm)
%   Acc:            Maximum distance between localised emitter from Spider
%                   and simulated emitter before it is considered a false
%                   positive
%   MapSize:        The size of the maps being used
%   Zoom:           The zoom being used for the SPIDER results in order to
%                   conver the results to nm
% 
% 
% Output:
% 
%   Recall:          The percentage of emitters that were correctly
%                   localised
%   Accuracy:       The accuracy range of the localisations. I.e. the 
%                   average distance between the correctly localised 
%                   emitters and the simulated emitters (e.g. all emitters 
%                   were correctly localised within a range of 30nm)
%   FalsePositive:  The number of false positive emitters in the Spider
%                   routine (a false positive emmiter is an emitter that 
%                   was localised in the Spider routine, but is not 
%                   actually there)
%   Sparsity:       The sparsity of the results, shown in %
% 
% *************************************************************************
% *************************************************************************

% Conversion of the simulated emitters from pixels to nm.
spider(:,1) = spider(:,1)*PixelToNm / Zoom;
spider(:,2) = spider(:,2)*PixelToNm / Zoom;
position = position*PixelToNm;

% Getting the size of the original image.
n = size(position,1);
o = size(spider,1);

% Pre-allocation of matrices to improve speed of calculations.
m = zeros(n,1);
loc = zeros(n,1);
XPos = zeros(n,1);
YPos = zeros(n,1);
Number = zeros(n,1);
Weight = zeros(n,1);
TotalDistance = zeros(n,1);
Distance = zeros(o,1);
Pos = zeros(o,1);
FPMatrix = [];

% Setting values to 0.
FalsePositive = 0;
LocalisationAccuracy = 0;

% Calculation of the mean centered emitters. I.e. when more than one
% emitter has been localised for a given simulated emitter, it will
% calculate the 'average' emitter first. The localised emitters are 
% assigned to their closest simulated emitter.
% The emitter has to be localised within a certain radius of the simulated
% emitter before it is used for further calculations.
for i = 1:o
    
    for j = 1:n
        
        % Calculation of the distance between all the localised emitters
        % and the simulated emitters.
        m(j) = sqrt((spider(i,1) - position(j,1))^2 + (spider(i,2) - position(j,2))^2);
        
    end
    
    % The smallest distance is calculated and the localised emitter is
    % assigned to the corresponding simulated emitter.
    [Min,I] = min(m);
    
    % Check if this distance is smaller than the wanted range and assign to
    % this simulated emitter if this is the case. If it is not the case,
    % assign as a false positive emitter.
    if Min < Acc
        
        XPos(I) = XPos(I) + spider(i,2)*spider(i,3);
        YPos(I) = YPos(I) + spider(i,1)*spider(i,3);
        Number(I) = Number(I) + 1;
        Weight(I) = Weight(I) + spider(i,3);
        Pos(i) = I;
        
    else
        
        FPMatrix = [FPMatrix;spider(i,1) spider(i,2) spider(i,3)];
        FalsePositive = FalsePositive + 1;
        
    end
    
end

% Calculate the mean centered emitters.
MeanCentered = [(YPos ./ Weight) (XPos ./ Weight)];
MeanCentered = MeanCentered(isfinite(MeanCentered(:,1)),:);

% Calculate the distance between the localised mean centered emitters and
% the simulated emitters again to check the other statistics.
for i = 1:size(MeanCentered,1)
    
    for j = 1:n
        
        m(j) = sqrt((MeanCentered(i,1) - position(j,1))^2 + (MeanCentered(i,2) - position(j,2))^2);
        
    end
    
    [Min,I] = min(m);
    
    loc(I) = loc(I) + 1;
    LocalisationAccuracy = LocalisationAccuracy + Min;
        
end

% Calculate the emitters that are correctly localised.
% Calculate the localisation accuracy.
Recall = size(find(loc),1) / n * 100;
Accuracy = LocalisationAccuracy / size(find(loc),1);

% Calculate the sparsity of the results in %
Sparsity = o/n*100;

end