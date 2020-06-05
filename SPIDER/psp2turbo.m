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
% Very fast fitting of tensor P-splines to a full data matrix Z.
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

function F = psp2turbo(Z, Ppars, W);

% *************************************************************************
% *************************************************************************
% Input: 
%
%   Z:      Data matrix
%   Ppars:  Parameters of the fitting
%   W:  	Weights for the fitting
%
% Output:
%
%   F:      The fit of the data matrix Z
%
% *************************************************************************
% *************************************************************************

[m n] = size(Z);
if nargin < 3
    W = ones(m, n);
end

% Prepare bases and penalties
cpars = [0 n+1 Ppars(2, :)];
rpars = [0 m+1 Ppars(1, :)];
Bc = bsplbase((1:n)', cpars);
nc = size(Bc, 2);
Dc = diff(eye(nc), cpars(6));
Pc = cpars(5) * Dc' * Dc;
Br = bsplbase((1:m)', rpars);
nr = size(Br, 2);
Dr = diff(eye(nr), rpars(6));
Pr = rpars(5) * Dr' * Dr;

T = rowtens(Br, Br)' * W * rowtens(Bc, Bc);
T = reshape(T, [nr nr nc nc]);
T = permute(T, [1 3 2 4]);
T = reshape(T, nr * nc, nr * nc) + kron(Pc, eye(nr)) + kron(eye(nc), Pr);
s = reshape(Br' * (W .* Z) * Bc, nr * nc, 1);
A = reshape(T \ s,  nr, nc);
F = Br * A * Bc';
