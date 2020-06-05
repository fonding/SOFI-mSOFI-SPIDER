%stack=sofiDetect(stack,scale,noise,bias,bits)
%---------------------------------------------
%
%Convert accumulated/detected photon counts to detector values.
%
%Inputs:
% stack  Detected image energy [photons]
% scale  Conversion scale [units/photon]        {1}
% noise  Read noise standard deviation [units]  {1}
% bias   Bias/offset value [units]              {2*noise}
% bits   Number of bits in converted values     {16}
%
%Output:
% stack  Detector/sensor values

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
function stack=sofiDetect(stack,scale,noise,bias,bits)
if nargin < 2 || isempty(scale)
   scale=1;
end
if nargin < 3 || isempty(noise)
   noise=1;
end
if nargin < 4 || isempty(bias)
   bias=2*abs(noise);
end
if nargin < 5 || isempty(bits)
   bits=16;
end
if bias
   stack=single(stack) + bias;
end
if noise
   stack=single(stack) + noise*randn(size(stack),'single');
end
stack=max(0,min(round(stack),pow2(1,bits)-1));
if bits < 9
   stack=uint8(stack);
elseif bits < 17
   stack=uint16(stack);
end
