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
% 'Interface' that manages all the steps in the different calculations
% needed to obtain the sparse result of the data matrix.
% It is a two-step algorithm for the deconvolution of the different
% patches.
%
% Authors:  
%       - HUGELIER, S.      (1)
%       - DE ROOI, J.J.     (2)
%       - DEVOS, O.         (1)
%       - EILERS P.H.C.     (3)
%       - RUCKEBUSCH, C.    (1)
%           
% (1):  LAboratoire de Spectrochimie Infrarouge et Raman (LASIR)
%       Université de Lille 1, UMR CNRS 8516
%       Bât C5, Cité Scientifique
%       59655 Villeneuve d'Ascq - France%           
% (2):  Biosystems Data Analysis
%       Swammerdam Institute for Life Sciences (Universiteit van Amsterdam)
%       Room C2.205, Science Park 904
%       1098 XH Amsterdam - The Netherlands
% (3):  Department of Biostatistics
%       Erasmus Medical Center
%       Room Na-2418, Wytemaweg 80
%       3015 CN Rotterdam - The Netherlands
%
% *************************************************************************

function Xsave = SpiderCalculations(image,kappa, B, n2, S, B0, n0, n02, S0, sel, h, R0, R)

% *************************************************************************
% *************************************************************************
% Input: 
%
%   image:  Patch used for the calculation of the emitters
%   kappa:  The penaltyfactor that determines the sparsity of the image
%   
%   The other input variables are variables needed for the calculation but
%   need no further explanation as they are automatically calculated by the
%   SPIDER interface program.
%
% Output:
%
%   Xsave:  The sparse result of the patch after the calculations 
%
% *************************************************************************
% *************************************************************************

% Important parameters
% Weight of the penalty for zoom
kappa0 = 0.01;
kappa_pos = 1000000;
% Safety parameter
beta = 0.00001;

% Do multiple images
Xsave = zeros(n2, 1);
Y = image;
rhs0 = reshape(B0' * Y * B0, n02, 1);
    
% Solve without zoom
xhat0 = solve_L0pen(S0, R0, rhs0, kappa0, kappa_pos, beta, sel);
Xhat0 = reshape(xhat0, n0, n0);
Xhat = Xhat0(h, h);
xhat = Xhat(:);
sel = find(xhat > 0);
    
rhs = reshape(B' * Y * B, n2, 1);
xhat = zeros(n2, 1);

xhat = solve_L0pen(S, R, rhs, kappa, kappa_pos, beta, sel);
    
Xsave(:, 1) = xhat;
end

