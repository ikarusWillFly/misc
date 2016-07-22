function screen2png(filename, siz, h)
%SCREEN2PNG Generate a PNG file of the current figure with
% dimensions consistent with the figure's screen dimensions.
%
% SCREEN2PNG('filename') saves the current figure to the
% PNG file "filename".
%
% Sean P. McCarthy
% Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved

if nargin < 1
	error('Not enough input arguments!')
end
if nargin < 2,
	% window scaling factor of 1
	siz = 1;
end
if nargin < 3,
	% figure handle
	h = gcf;
end
if size(siz,2) == 3,
	dpi = siz(3);
	siz = siz(1:2);
else
	dpi = 100;
end

oldscreenunits = get( h, 'Units' );
oldpaperunits  = get( h, 'PaperUnits' );
oldpaperpos    = get( h, 'PaperPosition' );

set( h, 'Units', 'pixels' );
oldsiz = get( h, 'Position' );
switch numel(siz),
	case 1 % scale
		siz = oldsiz * siz;
	case 2 % size
		% checks if any of the inputs is NaN; calculate if needed
		switch sum([ 0 find(isnan(siz)) ]),
			case 0 % none
				% Nothing to do
			case 1
				siz(1) = siz(2) * oldsiz(3)/oldsiz(4);
			case 2
				siz(2) = siz(1) / (oldsiz(3)/oldsiz(4));
			otherwise
				error('Either width or height should be NaN, not both');
		end
		siz = [ 0 0 siz ];
	otherwise
		error('Size must be either [ X Y ] or a scaling factor');
end

set( h, 'PaperUnits', 'inches', 'PaperPosition', siz/100 )
print( h,'-dpng', filename, '-r100' );
drawnow
set( h, 'Units', oldscreenunits, 'PaperUnits', oldpaperunits, ...
	'PaperPosition', oldpaperpos )
