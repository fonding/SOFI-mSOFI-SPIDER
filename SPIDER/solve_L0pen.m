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
% Function used to fit the sparse data with the single emitters. It uses a
% very accurate approximation of the L0-norm by using an L2-norm. The
% values under a certain threshold value are removed (as they are no real
% emitters and will never be in subsequent iterations).
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
%       59655 Villeneuve d'Ascq - France
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

function xhat = solve_L0pen(S, R, rhs, kappa, kappa_pos, beta, sel)

% *************************************************************************
% *************************************************************************
% Input: 
%
%   S:          The inproduct of kron(B) with itself
%   R:          Sparse matrix that contains the penaltyfactor for every
%               pixel of the matrix (the same for every one)
%   rhs:        Reconvolved data matrix
%   kappa:      The penaltyfactor that determines the sparsity of the image
%   kappa_pos:  The penaltyfactor to exclude negative emitters from the
%               result
%   beta:       A safety factor in case the values in xhat get close to 0
%   sel:        Selection matrix that tells the program which values should
%               be refined and which ones should be excluded
%
% Output:
%
%   xhat:      The sparse matrix
%
% *************************************************************************
% *************************************************************************

 n = size(S, 2);
 xhat = zeros(n, 1);
    for it = 1:50
        xold = sparse(xhat);
        xh = (S(sel, sel) + R(sel, sel)) \ rhs(sel);
        xhat = sparse(0 * xhat);
        xhat(sel) = xh;
        r = sparse(1 ./ (xhat .^ 2 + beta ^ 2));
        R = sparse(diag(kappa * r + kappa_pos * (xhat<0)));
        dx = max(abs(xhat -xold));
        sel = find(r < 1000);
        if dx < 1e-4
            break
        end
    end
end