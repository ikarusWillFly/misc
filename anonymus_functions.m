
%%% functions
% fun allows you to perform an anonymous function over a variable. This
% might be useful in case of long variable names or structures:
fun         = @(fun,dat) cell2mat(cellfun(fun,{dat},'uniformoutput',false));
fun         = @(fun,dat) cell2mat(cellfun(fun,dat{:},'uniformoutput',false)); % this should work with variable number of inputs

%%% cell functions
% uncell extract cells on a specified dimension
uncell      = @(x,d) cat(d,x{:});
% uncellSz extract cells keeping the original dimensions of the data, NB:
% if the cell is a vector, matlab automatically adds a dimension
uncellSz    = @(x) reshape(cat(ndims(x{1}),x{:}),[size(x{1}),size(x)]);

%%% size functions
sizes       = @(data,dims) cellfun(@(x) size(data,x),num2cell(dims)); % like size, but for multiple dimensions
matchSize   = @(x,y) repmat(x,size(y)./size(x)); % it should replicate a matrix y to match the size of the matrix x
cutDim{1}   = @(x,y) x(y,:,:,:,:); cutDim{2} = @(x,y) x(:,y,:,:,:); cutDim{3} = @(x,y) x(:,:,y,:,:); % cutDim selects specific ranges of a matrix
cutCell     = @(data,dim) cellfun(@(cell) cutDim{dim}(cell,1:min(cellfun(size(data,dim)))),data,'uniformoutput',false); % it cuts any cell to the size of the smallest cell

find_col    = @(data) (find(data) - (rem(find(data)-1, size(data,1)) + 1))/size(data,1) + 1;
find_row    = @(data)  rem(find(data)-1, size(data,1)) + 1;

%%% transformations
roundTo     = @(number,precision) round(number*(1/precision))*precision; % rounds a number for a certain precision
norm01      = @(dat,dim) bsxfun(@rdivide,bsxfun(@minus,dat,min(dat,[],dim)),range(dat,dim)); % normalize between 0 and 1 for a certain dimention
closest     = @(x,y) mod(find(ismember(abs(bsxfun(@minus,x,y(:))),min(abs(bsxfun(@minus,x,y(:)))'))')-1,size(x,2))+1; % finds the index for the entry in a vector closest to a certain value 
mean_bin    = @(data,bins) cellfun(@(val)  mean(data(ismember(bins,val))), num2cell(unique(bins)));   % it averages the entries of the matrix data indexed by the entries of a second matrix bins
bin_fun     = @(fun,data,bins) cellfun(@(val) fun(data(ismember(bins,val))), num2cell(unique(bins))); % performs any function on the entries of the matrix data indexed by the entries of a second matrix bins
dimfun      = @(fun,dat,dim) uncellSz(cellfun(fun,num2cell(dat,dim),'uniformoutput',false)); % perform any function on a specific dimension of the data, good also for functions that don´t allow 2 or three dimension matrixes
isinside    = @(val,bnd) val > bnd(1) & val < bnd(2); % is the value val within the values (bnd)
isoutside   = @(x,y) x < y(1) | x > y(2); % is the value val outside the values (bnd)
isodd       = @(x) round(x./2) == (x./2); % is odd
logBase     = @(x,base) log10(x)./log10(base); % it calculates a logarithm of any base
unit        = @(x) 1*exp(1i*angle(x)); % transforms a complex vector in a unity vector
complex2real = @(x) abs(x).*cos(angle(x)); % trasforms a complex vecotr in real signal
rad2deg     = @(rad) (180/pi).*rad
%%% plotting
subIndPlot  = @(Yfig,Xfig,Yplot,Xplot)  subplot(Yfig,Xfig,Xplot + (Yplot-1)*Xfig); % finally you can specify the subplot based on the row and columns :D
fullFigure  = @() figure('position',get(0,'ScreenSize')); % creates a full screen figure
dockFigure  = @() figure('windowstyle','docked'); % plots a docked figure 
closeExcept = @(fig) close(setxor(findall(0,'type','figure'),fig)); % close all figures but the one specified
lineH       = @(y,spec) plot(repmat(xlim,numel(y),1)',repmat(y(:),1,2)',spec{:}); % plot 2d horizontal lines  IE: lineV([2:3:6],{':k','linewidth',3})
lineV       = @(x,spec) plot(repmat(x(:),1,2)',repmat(ylim,numel(x),1)',spec{:}); % plot 2d vertical lines    IE: lineV([2:3:6],{':k','linewidth',3})
lineZ       = @(x,y,spec) plot3(repmat(x(:),1,2)',repmat(y(:),1,2)',repmat(zlim,numel(x),1)',spec{:}); % plots 3D lines on spanning z axis in position x and y, IE: lineZ([2:3:6],0,0,{':k','linewidth',3})
lineY       = @(x,z,spec) plot3(repmat(x(:),1,2)',repmat(ylim,numel(x),1)',repmat(z(:),1,2)',spec{:}); % plots 3D lines on spanning y axis in position x and z, IE: lineY([2:3:6],0,0,{':k','linewidth',3})
lineX       = @(y,z,spec) plot3(repmat(xlim,numel(y),1)',repmat(y(:),1,2)',repmat(z(:),1,2)',spec{:}); % plots 3D lines on spanning x axis in position y and z, IE: lineX([2:3:6],0,0,{':k','linewidth',3})
centClim    = @(x) set(x,'clim',max(abs(get(x,'clim')))*[-1 1]); % centers the colorbar at 0 based on the absolute maximum of the data
toPatch     = @(x1,x2,y1,y2,spec) patch([x1(:);flipud(repmat(x2(:),numel(x1)./numel(x2),1))],[y1(:);flipud(repmat(y2(:),numel(y1)./numel(y2),1))],spec{:}); % create patches: IE: patch(0:10,0:10,ones(10,1)-1,1,{'r','facealpha',0.2})
patchCircle = @(diam,X,Y,spec) fun(@(p) patch(sin(p)*diam+X,cos(p)*diam+Y,spec{:}),-pi:.1:pi); % patches a circle of diameter DIAM, centered in x,y.
findaxes    = @(cf) findobj(allchild(cf),'type','axes');
setFields   = @(strct,field,value) arrayfun(@(strct) setfield(strct,field,value),strct,'UniformOutput', false);

imagesc_patch = @(data,alpha,spec) arrayfun(@(x,y,z) patch([0,1,1,0]-.5+y,[0,0,1,1]-.5+x,spec{:},'facealpha',min(z,1)),...
    rem(find(data)-1, size(data,1)) + 1,...
    (find(data) - (rem(find(data)-1, size(data,1)) + 1))/size(data,1) + 1,...
    repmat(alpha,numel(find(data))./size(alpha,1),1));

%%% analysis / fieldtrip 
build_time         = @(d) cellfun(@(x) (0:length(x)-1)/d.fsample,d.trial,'uniformoutput',false); % it creates a fake time based on the fieltrip data structure
build_sampleinfo   = @(d,cont)   bsxfun(@plus,[[0;cumsum(cellfun(@length,reshape(d.trial(2:end),[],1)))], cumsum(cellfun(@length,d.trial(:)))-1],(2-cont*2)*(1:numel(d.trial))'); % it creates the sampleinfo of the data, cont determines whether the data is continuous or not
chanType           = @(data,type,hdr) find(ismember(hdr.chantype(ismember(data.label,hdr.label)),type));  % find the channels in the data that belong to a certain type (EEG/EMG) in the hdr
log_foi            = @(fr_beg,fr_end,n_fr) logspace(log10(fr_beg),log10(fr_end),n_fr); % it creates a logarithmic spacing of frequencies for the analyisis
str_index          = @(old_index,new_index) cell2mat(cellfun(@(x) find(strcmpi(new_index,x)),old_index,'uniformoutput',false)); 
build_anova_factor = @(dat,dim) reshape(bsxfun(@times,ones(size(dat)),shiftdim([1:size(dat,dim)]',1-dim)),[],1); % it creates a dummy variable for the grouping of anova data (dat), given the structure of dat 1 dimension X 1 factor 
find_labels        = label(cellfun(@(x) ~isempty(strfind(x,'C'))&&~isempty(str2num(x(end)))&&mod(str2num(x(end)),2)==1,label));
build_freq_names   = @(foi,smoothing) strsplit(sprintf('%d-%d Hz,',cell2mat(arrayfun(@(x,y) x + y./2*[-1 1]',foi(:),repmat(smoothing(:),numel(foi)/numel(smoothing),1),'uniformoutput',false))),',');

%%% inline conditionals
inline_try = @(dat,fun,fun2)      evalc('try     fun(dat), catch fun2(dat), end;');
inline_if  = @(dat,cond,fun,fun2) evalc('if cond fun(dat), else  fun2(dat), end;');

%%% MATH
gauss       = @(param,xval) param(1)*exp(-((xval-param(2)).^2)/(2*param(3))^2);  % 1 = heigh; 2 = peak location; 3 = width;
sigm        = @(param,xval) 0+(param(1)-0)./(1+exp((param(2)-xval)*param(3))) ;  % 1 = heigh; 2 = threshold; 3 = slope;


%% ROUTINES

%%% open figure
close all; fullFigure(); ca = gca; map = ca.ColorOrder;

%%% set FFT analysis bands
bands              = [8 25; 15 25; 9 35;25 40];
[foi, ord]         = sort(mean(bands'));
smoothing          = range(bands(ord,:)');

%%% determine rows and columns subplots
n_subplots = factor(n_plots);
while numel(n_subplots) > 2
    n_subplots = sort([prod(n_subplots(1:2)),n_subplots(3:end)]);
end

%%  TRASH

%{
if false
    subjCode = @(x) x(1:strfind(x,'_')-1)
%% MISC
doNothing  = @(x) x;
fun        = @(fun,dat) cell2mat(cellfun(fun,{dat},'uniformoutput',false));

uncell     = @(x,d) cat(d,x{:});

norm01     = @(dat) reshape((dat-min(dat(:)))./range(dat(:)),size(dat));                    % normalize between 0 and 1 across all the matrix
norm01     = @(dat,dim) bsxfun(@rdivide,bsxfun(@minus,dat,min(dat,[],dim)),range(dat,dim)); % normalize between 0 and 1 for a certain dimention

matchSize  = @(x,y) repmat(x,size(y)./size(x));

%%% allows you to use cellfun without writing the uniformoutput input
dimfun     = @(dat,fun,dim) uncellSz(cellfun(fun,num2cell(dat,dim),'uniformoutput',false));
cellfun2   = @(x,y) cellfun(x,y,'uniformoutput',false);
cellfun2   = @(x,y) inline_try(x,@(x)cellfun(x,y),@(x)cellfun(x,y,'uniformoutput',false));

%%% allow for in line data cutting
cutDims    = @(x,y,z,j) x(y,z,j);
cutDim1    = @(x,y) x(y,:,:);     cutDim2   = @(x,y) x(:,y,:);     cutDim3   = @(x,y) x(:,:,y);
cutDim(1)  = @(x,y) x(y,:,:,:,:); cutDim(2) = @(x,y) x(:,y,:,:,:); cutDim(3) = @(x,y) x(:,:,y,:,:);
cutDim{1}  = @(x,y) x(y,:,:,:,:); cutDim{2} = @(x,y) x(:,y,:,:,:); cutDim{3} = @(x,y) x(:,:,y,:,:);

%%% smooths bidimensional matrixed mantaining the original size
smoothR    = @(x,y) reshape(smooth(x,y),size(x));
smooth2d   = @(x,y,y2) reshape(smooth(reshape(smooth(x,y),size(x))',y2),size(x'))';

%%% rounds the number to any precision
roundTo    = @(number,precision) round(number*(1/precision))*precision;
closest    = @(x,y) mod(find(ismember(abs(bsxfun(@minus,x,y)),min(abs(bsxfun(@minus,x,y))'))'),size(x,2));

%%% inline conditionals (sweet)
inline_try = @(dat,fun,fun2)      evalc('try     fun(dat), catch fun2(dat), end;');
inline_if  = @(dat,cond,fun,fun2) evalc('if cond fun(dat), else  fun2(dat), end;');

calc_CI     = @(dof,alpha) atanh(1-alpha.^(1./(fix(dof)-1)));

gauss       = @(param,xval) param(1)*exp(-((xval-param(2)).^2)/(2*param(3))^2);  % 1 = heigh; 2 = peak location; 3 = width;
sigm        = @(param,xval) 0+(param(1)-0)./(1+exp((param(2)-xval)*param(3))) ;  % 1 = heigh; 2 = threshold; 3 = slope;
%% PLOT 
centClim     = @(x) set(x,'clim',max(abs(get(x,'clim')))*[-1 1]);
subIndPlot    = @(Yfig,Xfig,Yplot,Xplot)  subplot(Yfig,Xfig,Xplot + (Yplot-1)*Xfig);

%%% set the axis of the plot to fit the data
setAxis      = @(x) eval('axis tight; set(gca,''ylim'',get(gca,''ylim'')+get(gca,''ylim'').*[.1 .1]);');
cLim         = @(clims) set(gca,'clim',clims);
lineO        = @(y,spec) plot(repmat(xlim,numel(y),1)',repmat(y(:),1,2)',spec{:});
lineV        = @(x,spec) plot(repmat(x(:),1,2)',repmat(ylim,numel(x),1)',spec{:});
toPatch      = @(x1,x2,y1,y2,spec) patch([x1(:);flipud(repmat(x2(:),numel(x1)./numel(x2),1))],[y1(:);flipud(repmat(y2(:),numel(y1)./numel(y2),1))],spec{:}); 

tick         = @(axis,n) set(gca,[axis,'tick'],min(get(gca,[axis,'tick'])):range(get(gca,[axis,'tick']))/n:max(get(gca,[axis,'tick'])));
tick         = @(axis,n,rounding) set(gca,[axis,'tick'],round((min(get(gca,[axis,'tick'])):range(get(gca,[axis,'tick']))/n:max(get(gca,[axis,'tick'])))*(1/rounding))*rounding);
tick         = @(axis,n,rounding) set(gca,[axis,'tick'],roundTo(min(get(gca,[axis,'tick'])):range(get(gca,[axis,'tick']))/n:max(get(gca,[axis,'tick']))),rouding);

fullFigure  = @() figure('position',get(0,'ScreenSize'));
dockFigure  = @() figure('windowstyle','docked');
min(max(findobj('Type','figure')+1),n_fig)

%% ANALYSIS
fake_sampleinfo = @(d,cont)   bsxfun(@plus,[[0;cumsum(cellfun(@length,reshape(d.trial(2:end),[],1)))], cumsum(cellfun(@length,d.trial(:)))-1],(2-cont*2)*(1:numel(d.trial))');
fake_time       = @(d) cellfun(@(x) (0:length(x)-1)/d.fsample,d.trial,'uniformoutput',false);

%%% cat cells at the length of the smallest cell (useful prior concatenating data cells) 
cutCell   = @(data) cellfun(@(cell) cell(1:min(cellfun(@length,data)),:),data,'uniformoutput',false);
cutCell   = @(data,dim) cellfun(@(cell) cutDim{dim}(cell,1:min(cellfun(@(x) size(x,dim),data))),data,'uniformoutput',false);

%%% append fields contained in different structures embedded in cells 
appField  = @(struct,field,dim) cell2mat(shiftdim(cellfun(@(x) x.(field),data,'uniformoutput',false),dim));

%%% reshape frequency data as tapersXtrialsXchanXfreqXtime
getTapers = @(data) reshape(data.fourierspctrm,[data.cumtapcnt(1),numel(data.time),size(data,2),size(data,3),size(data,4)]);
unit      = @(x) 1*exp(1i*angle(x));
getTapers = @(dataFFT) reshape(dataFFT.fourierspctrm,[dataFFT.cumtapcnt(1),size(dataFFT.cumtapcnt,1),size(dataFFT.fourierspctrm,2),size(dataFFT.fourierspctrm,3),size(dataFFT.fourierspctrm,4)]);
%% MEASUREMENTS

%%% execute command in the background
clear b c
job = batch('[b,c] = find(max(rand(a)));');
wait(job)
load(job)

%%% execute command in the background separating job and task creation (faster)
job = createJob('job');
fun = @() find(max(rand(a)));
t = createTask(job,@(x) find(max(rand(x))),2,{10});
submit(job);
wait(job);
out = job.getAllOutputArguments;

%%% set timer for actions
t1 = timer('TimerFcn',@(x,y) disp(toc),'StartDelay',2);
t2 = timer('TimerFcn',@(x,y) disp(toc),'StartDelay',4);
start(t1)
start(t2)
tic

%%% take the cell relative to the current condition in case there are more
%%% than one NB: different conditions are in different columns
setCondition   = @(field,trial) cell2var(field(:,min(trial,size(field,2))));

%%% input checking
field = 'field';        value = 'value';
if ~isfield(cfg,field), cfg.(field) = value; end

%%% exit if esc is pressed
[ifPress,~,key] = KbCheck(0); if ifPress,  key = lower(KbName(key)); switch key case('esc'); keep_looping = 0; end, end

%%% initialize BVR
BVR.obj           = actxserver('VisionRecorder.Application');
BVR.view          = @() BVR.obj.Acquisition.ViewData;
BVR.stop          = @() BVR.obj.Acquisition.StopRecording;
BVR.DC            = @() BVR.obj.Acquisition.DCCorrection;

%%% initialize TCP sender
TCP.host       = '134.2.117.140'; % IPv4 of current machine
TCP.port       = 3000;
TCP.obj        = tcpip(TCP.host, TCP.port, 'NetworkRole', 'Client');
TCP.init_stim  = @() fwrite(TCP.obj,'init');
TCP.start_stim = @() fwrite(TCP.obj,'start');
TCP.stop_stim  = @() fwrite(TCP.obj,'stop');
fopen(TCP.obj);

%%% initialize TCP reciver
TCP.host       = '0.0.0.0'; % accept from all IPs
TCP.port       = 3000;      % conection port
TCP.timeout    = .01;       % time to wait to fill buffer (in seconds)
TCP.bufferSize = 1;         % size of buffer (in bytes)
TCP.obj        = tcpip(TCP.host, TCP.port, 'NetworkRole', 'Server');
set(TCP.obj,'Timeout',TCP.timeout);
set(TCP.obj, 'InputBufferSize',TCP.bufferSize)
fopen(TCP.obj);
read_data = @() char(fread(TCP.obj)');


%% cross correlation (to test)
maxLag       = 25;
cellfun2     = @(x,y) cellfun(x,y,'uniformoutput',false);
cutDim1      = @(x,y,z) x(y:z,:,:); cutDim2 = @(x,y,z) x(:,y:z,:); cutDim3 = @(x,y,z) x(:,:,y:z);
tmp          = cellfun(@squeeze,num2cell(cat(3,dataTmp.trial{:}),[2 3]),'uniformoutput',false);
[val,ind]    = cellfun2(@(x) max(cutDim2(xcorr(x,maxLag,'coeff'),2,size(x,2)+1)),tmp);
ind          = cat(1,ind{sel_chans});
ind          = mode(ind)-maxLag;

tmp  = cellfun(@(x,y) x(sel_chans,round(size(x,2)/2 +[-size(x,2)/2+maxLag+1:size(x,2)/2-maxLag-1]-y)),dataTmp.trial,num2cell(ind),'uniformoutput',false);
tmp2 = squeeze(cutDim1(cat(3,tmp{:}),25,25));
plot(tmp2);





%% PLV

UNCELL   = @(x,d) cat(d,x{:});
CELLTR   = @(x) cellfun(@(x) x',x,'uniformoutput',false);
CELLHIL  = @(x) cellfun(@hilbert,x,'uniformoutput',false);
UNIT     = @(x) cellfun(@(x) 1*exp(1i*angle(x)),x,'uniformoutput',false);
P_SHIFT  = @(x,ch) cellfun(@(x) bsxfun(@(x,y) x.*conj(y),x(:,ch),x(:,setxor(1:size(x,2),ch))),x,'uniformoutput',false);

HIL      = @(x) CELLHIL(CELLTR(x.trial));
ITC      = @(x) abs(mean(UNCELL(UNIT(HIL(x)),3),3));
PLV      = @(x) abs(mean(UNCELL(UNIT(P_SHIFT(HIL(x),numel(prm.ch))),3),3));



end

%}
