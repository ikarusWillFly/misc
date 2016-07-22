

%% MISC
doNothing  = @(x) x;

uncell     = @(x,d) cat(d,x{:});

norm01     = @(dat) (dat-min(dat))/range(dat);
%%% allows you to use cellfun without writing the uniformoutput input
cellfun2   = @(x,y) cellfun(x,y,'uniformoutput',false);
cellfun2   = @(x,y) inline_try(x,@(x)cellfun(x,y),@(x)cellfun(x,y,'uniformoutput',false));

%%% allow for in line data cutting
cutDims    = @(x,y,z,j) x(y,z,j);
cutDim1    = @(x,y) x(y,:,:);     cutDim2   = @(x,y) x(:,y,:);     cutDim3   = @(x,y) x(:,:,y);
cutDim(1)  = @(x,y) x(y,:,:,:,:); cutDim(2) = @(x,y) x(:,y,:,:,:); cutDim(3) = @(x,y) x(:,:,y,:,:);

%%% smooths bidimensional matrixed mantaining the original size
smooth2    = @(x,y) reshape(smooth(x,y),size(x));
smooth2D    = @(x,y,y2) reshape(smooth(reshape(smooth(x,y),size(x))',y2),size(x'))';

%%% rounds the number to any precision
roundTo    = @(number,precision) round(number*(1/precision))*precision;

%%% inline conditionals (sweet)
inline_try = @(dat,fun,fun2)      evalc('try     fun(dat), catch fun2(dat), end;');
inline_if  = @(dat,cond,fun,fun2) evalc('if cond fun(dat), else  fun2(dat), end;');


%% PLOT

%%% set the axis of the plot to fit the data
setAxis      = @(x) eval('axis tight; set(gca,''ylim'',get(gca,''ylim'').*[1.1 1.1]);');
cLim         = @(clims) set(gca,'clim',clims);
%% ANALYSIS
%%% cat cells at the length of the smallest cell (useful prior concatenating data cells) 
cutCell   = @(data) cellfun(@(cell) cell(1:min(cellfun(@length,data)),:),data,'uniformoutput',false);

%%% append fields contained in different structures embedded in cells 
appField  = @(struct,field,dim) cell2mat(shiftdim(cellfun(@(x) x.(field),data,'uniformoutput',false),dim));

%%% reshape frequency data as tapersXtrialsXchanXfreqXtime
getTapers = @(data) reshape(data.fourierspctrm,[data.cumtapcnt(1),numel(data.time),size(data,2),size(data,3),size(data,4)]);

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
ind = cat(1,ind{sel_chans});
ind = mode(ind)-maxLag;

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













