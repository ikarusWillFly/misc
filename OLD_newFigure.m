function handle = newFigure(handle)
if nargin < 1,
    figHandles = findobj('Type','figure');
    handle     = max([figHandles+1;1]);
end
figure(handle)
clf
set(handle,'windowstyle','docked')
figure(handle)
end

