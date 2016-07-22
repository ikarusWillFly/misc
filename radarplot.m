function ca = radarplot(cfg,data)
% ca = radarplot(cfg,data)
% INPUTS:
% data               = % data to plot matrix: Ncomponets X Nparameters 
% cfg.maxLims        = % MAXIMUM AXIS LIMITS       - double
% cfg.labels         = % LABELS OF THE PARAMETERS  - cell vector of chars: 1 X Ncomponents
% cfg.texts          = % TEXTS ON THE PATCHES EDGES - cell vector of chars: 1:Ncomponets X 1:Nparameters
% cfg.patchColors    = % COLORS OF THE PATCHES      - cell vector: Ncomponents X 1
% cfg.patchSpecs     = % PATCHES SPECIFICATIONS     - cell vector: Ncomponents X 1
% cfg.labelSpecs     = % LABELS SPECIFICATIONS      - cell vector: Ncomponents X 1
% cfg.textSpecs      = % EDGE TEXT SPECIFICATIONS   - cell vector: 1 X Nparameters
% 
% % DEFAULTS:
% defaultMaxLims     = max(abs(data(:)))*1.1;
% defaultLabels      = cellfun(@num2str,num2cell(1:n_prm),'uniformoutput',false);
% defaultTexts       = {};
% defaultPatchSpecs  = cellfun(@(x) {x,'faceAlpha',.2,'edgealpha',1,'edgeColor',x,'lineWidth',2},defaultPatchColors,'UniformOutput',false);
% defaultPatchColors = num2cell(ca.ColorOrder(1:n_cmp,:),2);
% defaultLabelSpecs  = {'k','fontSize',30};
% defaultTextSpecs   = num2cell(ca.ColorOrder(1:n_prm,:),2)';
% 
% % EXAMPLE OF USAGE
%
% clf;
% data             = rand(3,5); 
% cfg = [];
% cfg.maxLims      = max(data(:))*1.1;
% cfg.labels       = cellfun(@num2str,num2cell(1:size(data,2),1),'uniformoutput',false);
% cfg.texts        = repmat({'a','b','c','d','e'},size(data,1),1)
% cfg.patchColors  = {'r';'g';'b'};
% cfg.patchSpecs   = cellfun(@(x) {x,'faceAlpha',.2,'edgealpha',1,'edgeColor',x,'lineWidth',2},cfg.patchColors,'UniformOutput',false);
% cfg.labelSpecs   = {'k','fontSize',30};
% cfg.textSpecs    = {{'fontSize',20}};
% ca = radarplot(cfg,data)
% 


%%
ca    = gca;
n_cmp = size(data,1);
n_prm = size(data,2);

defaultMaxLims     = max(abs(data(:)))*1.1;
defaultLabels      = cellfun(@num2str,num2cell(1:n_prm),'uniformoutput',false);
defaultTexts       = {};
defaultPatchColors = num2cell(ca.ColorOrder(1:n_cmp,:),2);
defaultPatchSpecs  = cellfun(@(x) {x,'faceAlpha',.2,'edgealpha',1,'edgeColor',x,'lineWidth',2},defaultPatchColors,'UniformOutput',false);
defaultLabelSpecs  = {'k','fontSize',30};
defaultTextSpecs   = num2cell([repmat([{'fontSize'},{22},{'color'}],n_prm,1),num2cell(ca.ColorOrder(1:n_prm,:),2)],2);
defaultAxisSpecs   = {'w','edgecolor','k','linewidth',3};

field = 'labels';             value = defaultLabels;
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'texts';              value = defaultTexts;
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'patchSpecs';         value = defaultPatchSpecs;
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'patchColors';        value = defaultPatchColors;
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'textSpecs';          value = defaultTextSpecs;
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'maxLims';            value = defaultMaxLims;
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'labelSpecs';        value  = defaultLabelSpecs;
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'axisSpecs';        value   = defaultAxisSpecs;
if ~isfield(cfg,field), cfg.(field) = value; end

%% plot patches and texts
% the function "select_parameter" allows you to get a specific parameter
% for a certain condition in case it is specified. 
% In case of a common parameter for different conditions,
% the function returns the common value without returning matrix errors  
select_parameter = @(x,row,col) x{min(size(x,1),row), min(size(x,2),col)}; % get the parameter relative to the condition specified by row and column, in case it exist.

tmpDat        = [data data(:,1)];                                          % repeat the coordinates of the last point of the data to close the shape
tmpDat        = cat(1,repmat(cfg.maxLims,1,size(tmpDat,2)),tmpDat);        % treat the axis as another patch

% find the coordinates of the vertices and the texts
meshIn        = 2*pi/n_prm*[0:n_prm]+pi/n_prm;
[Theta, ~]    = meshgrid(meshIn,ones(1,size(tmpDat,1)));

X             = tmpDat.*sin(Theta);
Y             = tmpDat.*cos(Theta);
Xtext         = tmpDat*1.1.*sin(Theta);
Ytext         = tmpDat*1.1.*cos(Theta);

% concatenate the specifications for the components with the ones for the axis
texts         = cat(1,cfg.labels,cfg.texts);
textSpecs     = cat(1,{cfg.labelSpecs},cfg.textSpecs);
patchSpecs    = cat(1,{cfg.axisSpecs},cfg.patchSpecs);

i_cmp = 1; i_prm = 1;
for i_cmp     = 1 : n_cmp+1 % for the number of components (plus background)
    patchSpec = select_parameter(patchSpecs,i_cmp,1); % in case of a different patch specification for each component or parameter, select that one.
    patch(X(i_cmp,:),Y(i_cmp,:),patchSpec{:})
    
    if size(texts,1)>= i_cmp % if there is a specified text for that component
        for i_prm = 1 : n_prm
           textSpec = select_parameter(textSpecs,i_cmp,i_prm); % in case of a different text specification for each component or parameter, select that one.
           text(Xtext(i_cmp,i_prm),Ytext(i_cmp,i_prm),texts(i_cmp,i_prm),textSpec{:},'HorizontalAlignment','center');
        end
    end
end
axis equal; axis off
