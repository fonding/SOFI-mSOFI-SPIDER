% *************************************************************************
% SParse deconvolution of hIgh-Density supER-resolution images 
% *************************************************************************
%
% Reference to the publication:
%   Hugelier, S., de Rooi, J.J., Bernex, R., Duw?, S., Devos, O., Sliwa, M.,
%   Dedecker, P., Eilers, P.H.C. & Ruckebusch, C. (2016). Sparse 
%   deconvolution of high-density super-resolution images, Sci. Rep. 6, 
%   21413; doi: 10.1038/srep21413 (2016)
%
% Computes the (flattened) tensor products per row of the input.
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

function C = rowtens(A, B)

% *************************************************************************
% *************************************************************************
% Input: 
%
%   A & B:  Data matrices needed to do the fitting of the data matrix Z for
%           detection of the emitters
%
% Output:
%
%   C:      Tensor product of the input 
%
% *************************************************************************
% *************************************************************************

na = size(A, 2);
nb = size(B, 2);
ea = ones(1, na);
eb = ones(1, nb);
C = kron(A, eb) .* kron(ea, B);
