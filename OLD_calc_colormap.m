% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % cmap    = calc_colormap(tail, contrast, mode, resolution)
% % tail == 2, bothsided
% % tail == -1, blue
% % tail == 1, red
% % contrast = 0, no contrast
% % contrast >0 & <1, amount of backgrounded values
% % mode: background color: 'black,'white','gray','green'
% % resolution = fineness of the steps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap    = calc_colormap(tail, contrast, mode, resolution, saturation)

if nargin <1, tail = 2; end
if nargin <2, contrast = 0; end
if nargin <3, mode = 'black'; end
if nargin <4, resolution = 100; end
if nargin <5, 
    if strcmp(mode,'white') || strcmp(mode,'black') || strcmp(mode,'cont'),
        saturation = 0; 
    elseif strcmp(mode,'green') || strcmp(mode,'gray'),
        saturation = 0.5; 
    end
end

if contrast <0, contrast = 0; end
if contrast >1, contrast = 1; end

saturation  = sqrt(saturation);
resolution  = resolution./abs(tail);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construction of of elements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(mode,'black'),

grad        = linspace(1,0,resolution.*(1-contrast))';
null        = linspace(saturation,0,resolution.*(1-contrast))';
fill        = zeros(fix(resolution.*contrast),3);
blue        = cat(1,cat(2,null,null,grad),fill);
red         = cat(1,cat(2,grad,null,null),fill);

end

if strcmp(mode,'white'),

grad        = linspace(1,1,resolution.*(1-contrast))';
null        = linspace(0,(1-saturation),resolution.*(1-contrast))';
fill        = ones(fix(resolution.*contrast),3);
blue        = cat(1,cat(2,null,null,grad),fill);
red         = cat(1,cat(2,grad,null,null),fill);

end

if strcmp(mode,'green'),

grad        = linspace(1,saturation,resolution.*(1-contrast))';
null        = zeros(size(grad));
gren        = (cat(2,zeros(1,floor((resolution.*(1-contrast)./2))),linspace(0,saturation,ceil(resolution.*(1-contrast)./2)))');
if size(gren,1) < size(grad,1), gren = cat(1,0,gren); end
fill        = cat(2,zeros(floor(resolution.*contrast),1),ones(floor(resolution.*contrast),1).*saturation,zeros(floor(resolution.*contrast),1));
blue        = cat(1,cat(2,null,gren,grad),fill);
red         = cat(1,cat(2,grad,gren,null),fill);

end

if strcmp(mode,'gray'),
    
grad        = linspace(1,saturation,resolution.*(1-contrast))';
gray        = cat(2,zeros(1,floor((resolution.*(1-contrast)./2))),linspace(0,saturation,ceil(resolution.*(1-contrast)./2)))';
if size(gray,1) < size(grad,1), gray = cat(1,0,gray); end
fill        = ones(fix(resolution.*contrast),3).*saturation;
blue        = cat(1,cat(2,gray,gray,grad),fill);
red         = cat(1,cat(2,grad,gray,gray),fill);
end

if strcmp(mode,'cont'),
    
grad        = linspace(1,1,resolution.*(1-contrast))';
null        = linspace(0,(1-saturation),resolution.*(1-contrast))';
fill        = ones(fix(resolution.*contrast),3);
blue        = cat(1,cat(2,null,null,grad),fill);
red         = cat(1,cat(2,grad,null,null),fill);
red         = cat(1,red(1,:),red);
blue        = cat(1,blue(1,:),blue);
%red(:,2)    = linspace(0.5,0,length(red));
%blue(:,2)   = linspace(0.5,0,length(red));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construction of colormap    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if tail ==2,
    cmap = (cat(1,blue,flipud(red)));
end
if tail == -1,
    cmap = blue;
end
if tail == 1,
    cmap = red;
end