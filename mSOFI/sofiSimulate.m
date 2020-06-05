%stack=sofiSimulate(frames,centers,Ion,Ton,Toff,fwhm,Ibg)
%--------------------------------------------------------
%
%Simulate the acquisition of an image sequence of blinking emitters.
%
%Inputs:
% frames    Number of images
% centers   Emitter coordinates [x y] [pixel]
% Ion       Signal per frame (on state) [photon]
% Ton       Average duration of the on state [frame]
% fwhm      Point-spread function full width at half-maximum [pixel]
% Ibg       Background intensity [photon] {none}
%
%Output:
% stack     Image sequence [101 x 101 x frames]

%Copyright ? 2012 Marcel Leutenegger et al, École Polytechnique Fédérale de Lausanne,
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
function stack=sofiSimulate(frames,centers,Ion,Ton,Toff,fwhm,Ibg)
emitters=size(centers,1);
Ion=repmat(Ion(:),emitters/numel(Ion),1);
Ton=repmat(Ton(:),emitters/numel(Ton),1);
Toff=repmat(Toff(:),emitters/numel(Toff),1);
if nargin < 7 || isempty(Ibg)
   Ibg=0;
end
fig=statusbar('Photons...');
photons=zeros(emitters,frames);
for emitter=1:emitters
   photons(emitter,:)=Ion(emitter)*brightness(frames,Ton(emitter),Toff(emitter));
   fig=statusbar(emitter/emitters,fig);
   if isempty(fig)
      return;
   end
end
y=-50:50;
n=numel(y);
x=gaussian(fwhm,centers(:,1),y);
y=gaussian(fwhm,centers(:,2),y).';
fig=statusbar('Images...',fig);
stack=uint16(0);
stack(n,n,frames)=0;
n=ones(1,n,'uint8');
for frame=1:frames
   stack(:,:,frame)=poissrnd(x*(y.*photons(:,frame(n)))+Ibg);
   fig=statusbar(frame/frames,fig);
   if isempty(fig)
      stack=stack(:,:,1:frame);
      break;
   end
end
delete(fig);


%Noise-free image sequence of a blinking emitter.
%
% Nf        Number of frames
% Ton       On state lifetime [frames]
% Toff      Off state lifetime [frames]
% photons   Intensity trace [photons/frame]
%
function photons=brightness(Nf,Ton,Toff)
cycle=Ton + Toff;          % average blinking cycle time
cycles=10 + ceil(Nf/cycle);
times=[-Toff*log(rand(1,cycles));-Ton*log(rand(1,cycles))];
if rand < Ton/cycle        % initialize with a random start state:
   times(1)=0;             % exponential distribution has no memory
end
times=cumsum(times(:));    % times of on- and off-switching events
while times(end) < Nf
   cycles=ceil(2*(Nf - times(end))/cycle);
   cycles=[-Toff*log(rand(1,cycles));-Ton*log(rand(1,cycles))];
   cycles(1)=cycles(1) + times(end);
   times=[times;cumsum(cycles(:))];
end
times=times.';
Ton=times(2:2:end) - times(1:2:end);
photons=[zeros(size(Ton));Ton];
photons=cumsum(photons(:));
photons=diff(interp1(times,photons,0:Nf,'linear',0));


%One-dimensional Gaussian point spread function.
%
% fwhm      Full width at half-maximum diameters [pixel]
% centers   Emitter coordinates [pixel]
% points    Pixel coordinates
% psf       Gaussian PSFs
%
function psf=gaussian(fwhm,centers,points)
psf=fwhm/sqrt(8*log(2));
if numel(psf) > 1
   [~,psf]=ndgrid(points,psf);
end
points=points(:);
points=[[1.5 -0.5]*points(1:2);conv(points,[0.5;0.5],'valid');[-0.5 1.5]*points(end-1:end)];
[points,centers]=ndgrid(points,centers);

psf=diff(normcdf(points,centers,psf),[],1);
