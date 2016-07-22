% makes plots of potentials on head
% usage showfield(z,loc,scale, markersize,resolution, cbar) ;
%
% input:
%  z  values of field/potential to be shown.
%  loc   matrix containing 2D coordinates of channels in second and third column 
%  scale sets the scale of the color map. Either a 1x2 vector
%            corresponding to minimum and a maximum or just a number (say x) 
%            then the scale is from [-x x]. The default is 
%            {=[ -max(abs(z)) max(abs(z)) ]}
%  markersize is the size of the markers for the original channels
%  resolution sets the resolution, {100}
%  cbar=1 draws colorbar, otherwiseno colorbar is not shown. {1}
%   rewritten by R.Bauer & M.Vukelic for CIN AG NPT 22.06.2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function showfield(z,loc,scale,markersize,resolution, interpolmode, cbar) 

if nargin <3, scale=[-max(max(abs(z))),max(max(abs(z)))]; end;
if nargin <4, markersize = 24; end;
if nargin <5, resolution=100; end;
if nargin <6, interpolmode = 1; end;
if nargin <7, cbar=0; end;

if length(scale) == 1, scale = [min([scale, -scale]) max([scale, -scale])]; end

[~,m]=size(loc);
if m==2;
  x=loc(:,1);
  y=loc(:,2);
else
  x=loc(:,2);
  y=loc(:,3);
end

xlin = linspace(1.4*min(x),1.4*max(x),resolution);
ylin = linspace(1.4*min(y),1.4*max(y),resolution);
[X,Y] = meshgrid(xlin,ylin);

Z = griddata(x,y,z,X,Y,'nearest');
if interpolmode == 1,
    Z = smooth2D(Z,repmat(fix(resolution/5),1,2),2);   
end

  % Take data within head
  rmax=1.02*max(sqrt(x.^2+y.^2));
  mask = (sqrt(X.^2+Y.^2) <= rmax);
  Z(mask == 0) = NaN;
  
  
surface(X,Y,zeros(size(Z)),Z,'edgecolor','none');shading interp;

%caxis([ - max(abs(z)) max(abs(z))]);
caxis(scale);
hold on;
plot(x,y,'.k','markersize',markersize);


%meanx=mean(loc(:,2)*.85+.45);
%meany=mean(loc(:,3)*.85+.45);
scalx=1;
drawhead(0,.0,rmax,scalx);
set(gca,'visible','off');

axis([-1.2*rmax 1.2*rmax -1.0*rmax 1.4*rmax]);
axis equal;
%axis([-1.4*rmax 1.4*rmax -1.0*rmax 1.4*rmax]);
if cbar==1
  h=colorbar;set(h,'fontweight','bold')
end

%plot(.985*rmax*sin((0:1000)/1000*2*pi), .985*rmax*cos((0:1000)/1000*2*pi),'linewidth',2,'color','k'); 
return; 



function drawhead(x,y,size,scalx)

cirx=(x+scalx*size*cos((1:1000)*2*pi/1000) )';ciry=(y+size*sin((1:1000)*2*pi/1000))';

plot(cirx,ciry,'k','linewidth',1);
hold on;

ndiff=20;
plot( [x  cirx(250-ndiff) ],[y+1.1*size ciry(250-ndiff)],'k','linewidth',1);
plot( [x  cirx(250+ndiff) ],[y+1.1*size ciry(250+ndiff)],'k','linewidth',1);


return;
