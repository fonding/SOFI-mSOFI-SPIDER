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
% Compute the inproduct of the Kronicker product of B (the B-spline basis)
% with itself, by using an array trick
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

function S = fast_inprod(B)

% *************************************************************************
% *************************************************************************
% Input: 
%
%   B:      The B-spline basis
%
% Output:
%
%   S:      The inproduct of kron(B) with itself
%
% *************************************************************************
% *************************************************************************

[m n] = size(B);
Bb = rowtens(B, B);
W = ones(m, m);
T = Bb' * W * Bb;
T = reshape(T, [n n n n]);
T = permute(T, [1 3 2 4]);
S = reshape(T, n * n, n * n);
