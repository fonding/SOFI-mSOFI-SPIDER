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
% Function that performs a fast computation of a B-spline basis of degree
% 'deg' at positions 'x' on a uniform grid with 'ndx' intervals between
% 'x0' and 'x1'.
% The original function was designed by Paul H. C. Eilers (1996)
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

function B = bsplbase(x, bpars)

% *************************************************************************
% *************************************************************************
% Input: 
%
%   x:      The positions at which the B-splines should be created
%   bpars:  The parameters for the creation of the B-splines
%
% Output:
%
%   B:      The B-spline basis
%
% *************************************************************************
% *************************************************************************

x0 = bpars(1); x1 = bpars(2); ndx = bpars(3); deg = bpars(4);
x = x(:);
if (min(x) < x0) | (max(x) > x1)
  disp('Some elements of x out of bounds !!')
  return
end
dx = (x1 - x0) / ndx;
t = x0 + dx * ((-deg):(ndx-1));
T = ones(size(x)) * t;
X = x * ones(size(t));
D = (X - T) / dx;
B = diff([zeros(size(x)), D <= 1]')';
r = [2:length(t) 1];
for k = 1:deg
  B = (D .* B + (k + 1 - D) .* B(:, r)) / k;
end;


