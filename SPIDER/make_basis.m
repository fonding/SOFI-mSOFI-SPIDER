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
% Function to create a Gaussian convolution on a zoomed scale
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

function B = make_basis(m, zoom, sigma)

% *************************************************************************
% *************************************************************************
% Input: 
%
%   m:      Size of the basis (without zoom)
%   zoom:   Zoom for the basis
%   sigma:  the sigma of the PSF of the Gaussian function
%
% Output:
%
%   B:      The Gaussian convolution basis
%
% *************************************************************************
% *************************************************************************

n = m * zoom;
sigmaz = sigma * zoom;
uz = (-n):n;
vz = exp(-(uz' / sigmaz) .^ 2 / 2);

% Basis on zoomed scale
Bz = zeros(n, n);
nv = length(vz);
for i = 1:n
    l = max(i - n, 1);
    u = min(i + n, n);
    l2 = max(l - i + n+1, 1);
    u2 = min(u - i + n+1, nv);
    Bz(i, l : u) = vz(l2 : u2);
end
n2 = n ^ 2;

% Coarsened basis
B = zeros(m, n);
k = 0;
for i = 1:m
    for j = 1:zoom
        k = k + 1;
        B(i, :) = B(i, :) + Bz(k, :);
    end
end
B = B / zoom;

