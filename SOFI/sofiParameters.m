%[ratio,density,bright]=sofiParameters(sofi,cond,order,average)
%--------------------------------------------------------------
%
%Estimate emitter parameters from flat cumulants.
%
%Inputs:
% sofi      Flat cumulants of orders 2 to 4
% cond      Illumination and condition dependent exponent.  {2}
%           Zero in the single molecule limit; otherwise
%           1.5 for TIRF or 2 for full-field illumination.
% order     Order at which parameters are sampled.          {2}
% average   Weighted average parameters if set.             {true}
%
%Outputs:
% ratio     On-time ratio
% density   Emitter density
% bright    Emitter brightness

%Copyright © 2012 Marcel Leutenegger et al, École Polytechnique Fédérale de Lausanne,
%Laboratoire d'Optique Biomédicale, BM 5.142, Station 17, 1015 Lausanne, Switzerland.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [ratio,density,bright]=sofiParameters(sofi,cond,order,average)
if nargin < 2 || isempty(cond)
   cond=2;
end
if nargin < 3 || isempty(order)
   order=2;
end
if numel(sofi) < 4 || any(cellfun(@isempty,sofi(2:4)))
   error('sofi:orders','Require cumulant images of 2nd, 3rd and 4th order.');
end
%
% Resample cumulants to order.
%
if order > 1
   [x,y]=xygrid(size(sofi{order}));
   kappa2=interpolate(sofi{2}(:,:,:),x,y);
   kappa3=interpolate(sofi{3}(:,:,:),x,y);
   kappa4=interpolate(sofi{4}(:,:,:),x,y);
else
   for n=2:4
      [x,y,z]=size(sofi{n});
      kappa4=sofi{n}([1:x x-1:-1:x-n+1],[1:y y-1:-1:y-n+1],:);
      kappa4=reshape(kappa4,n,ceil(x/n),n,ceil(y/n),z);
      kappa4=squeeze(sum(sum(kappa4,1),3)/n^2);
      eval(sprintf('kappa%d=kappa4;',n));
   end
end
%
% Estimate the emitter parameters.
%
%  Cumulant model: kappan = E{PSF^n} Ns Ion^n fn(rho)
%
%  On-time ratio polynomials fn(rho):
%
%     f1 = rho
%     f2 = rho(1-rho)
%     f3 = rho(1-rho)(1-2rho)
%     f4 = rho(1-rho)(1-6rho+6rho²)
%
%  Cumulant ratios:
%
%          kappa3        cond
%     K1 = ------ = (3/2)     Ion (1-2rho)
%          kappa2
%
%          kappa4        cond    2            2
%     K2 = ------ = (4/2)     Ion (1-6rho+6rho )
%          kappa2
%
kappa3=kappa3./kappa2/(3/2)^cond;
kappa4=kappa4./kappa2/(4/2)^cond;
bright=sqrt(max(0,3*kappa3.^2 - 2*kappa4));
ratio=max(0,min(0.5-0.5*kappa3./bright,1));
density=max(0,kappa2./(bright.^2.*ratio.*(1-ratio)));
%
% Weighted average using the 2nd-order cumulant as weight.
%
if (nargin < 4 || average) && size(kappa2,3) > 1
   kappa2=abs(kappa2);
   kappa3=1./sum(kappa2,3);
   ratio=sum(kappa2.*ratio,3).*kappa3;
   density=sum(kappa2.*density,3).*kappa3;
   bright=sum(kappa2.*bright,3).*kappa3;
end
