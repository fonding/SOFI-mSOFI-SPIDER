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
% Function to remove the baseline of the data. It uses an asymmetric
% smoothing of the baseline surface by using 2D P-splines.
%
% Authors:  
%       - DE ROOI, J.J.     (1)
%       - EILERS P.H.C.     (2)
%           
% (1):  Biosystems Data Analysis
%       Swammerdam Institute for Life Sciences (Universiteit van Amsterdam)
%       Room C2.205, Science Park 904
%       1098 XH Amsterdam - The Netherlands
% (2):  Department of Biostatistics
%       Erasmus Medical Center
%       Room Na-2418, Wytemaweg 80
%       3015 CN Rotterdam - The Netherlands
%
% *************************************************************************

function [Z,W] = baselineSpider(Y, p, nseg, lambda, WeightBackground)

% *************************************************************************
% *************************************************************************
% Input: 
%
%   Y:      Data of which the baseline should be subtracted
%   p:      Asymmetry factor (close to 0 for removal of positive peaks and
%           close to 1 for removal of negative peaks)
%   nseg:   The number of splines used for the baseline subtraction
%   lambda: The penaltyfactor used for the smoothing of the baseline
%
% Output:
%
%   Z:      The baseline subtracted data
%
% *************************************************************************
% *************************************************************************

if nargin < 5
    W = 0 * Y + 0.5;
else
    W = WeightBackground;
end
Ppars = [nseg 3 lambda 3; nseg 3 lambda 3];

for it = 1:20
    Z = psp2turbo(Y, Ppars, W);
    R = Y - Z;
    Wo = W;
    W = p .* (R > 0 ) + (1 - p) .* (R < 0);
    dw = sum(W(:) ~= Wo(:));
    if dw == 0
        break
    end
end
