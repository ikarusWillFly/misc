


%%
minFr   = 4;
maxFr   = 50;
frStep  = 1;
numFr   = 40;
foi     = logspace(log10(minFr),log10(maxFr),numFr);
Qi      = 3;

ep_size = .05;
clear DAT
for sb = 1 : numel(subject)
    subj    = subject{sb};
    PTH.sbj = [PTH.data,subj,PTH.sep];
    tmp     = load([PTH.sbj,'dataSource'],'dataSource');
    data    = tmp.dataSource; clear tmp
    fs      = data.fsample;
    for cnd = 1 : 2
        sel_trials = find(data.trialinfo==cnd);
        clear tmp dataFilt dataHil dataSlide data0 dataUnwrap datOut
        for fr = 1 : numel(foi)
            clc, fprintf('starting frequency n° %d of %d\n',fr,numel(foi))
            cfg = [];
            cfg.bpfilter    = 'yes';
            cfg.bpfilttype  = 'fir';
            cfg.trials      = sel_trials;
            cfg.bpfreq      = foi(fr) + foi(fr)./([-1 1]*Qi)./2;
            tmp             = ft_preprocessing(cfg,data);
            cfg = []; 
            cfg.detrend     = 'yes';
            cfg.resamplefs  = 1/0.005; 
            tmp             = ft_resampledata(cfg,tmp);
            dataFilt(:,:,:,fr)        = cat(3,tmp.trial{:});
            
        end
        dataFilt(:,:,:,end+1)         = repmat(sin(linspace(0,20*range(data.time{1})*pi,size(dataFilt,2))),size(dataFilt).* [1,0,1,0] + [0 1 0 1]);
        dataFilt = permute(dataFilt,[2 3 1 4]); % time X trial X channel X frequency
        dataHil  = cell2mat(cellfun(@(x) angle(hilbert(x)), num2cell(dataFilt,[1 2]),'UniformOutput', false));
%         clear dataFilt
        %%%  CUT TRIALS IN SLIDING WINDOWS
        cfg = [];
        cfg.time_dimension  = 1;                 % time dimension
        cfg.trial_dimension = 2;                 % trial dimension
        cfg.epoch_size      = round(ep_size*fs); % desired epoch sizes (in samples)
        cfg.sliding         = round(ep_size*fs); % desired sliding size (in samples)
        cfg.keeptrials      = false;             % flag to keep the trials separated
        cfg.feedback        = true;              % flag to receive feedback with datadimensions
        dataSlide           = slidingEpochs(cfg,dataHil);
        dataUnwrap          = unwrap(dataSlide);
%         clear dataSlide
        %%% MANIPULATE THE DATA
        data0               = bsxfun(@minus,dataUnwrap,dataUnwrap(1,:,:,:)); % baseline the epochs to the first phase
        datOut              = permute(data0,[2 1 3 4]);                      % get the repetitions dimension first
        datOut              = squeeze(std(datOut));                          % perform the standard deviation
        datOut              = datOut/(2*pi);                                 % calculate the percentage of a cycle equal to 2 std of the phases
%         clear data0 dataUnwrap
%         DAT(sb).(['tr_',num2str(cnd)]) = datOut;
        DAT{cnd} = datOut;
    end
    save([PTH.sbj,'phaseStab_Qi',num2str(Qi)],'DAT')
end


%%
size(dataSlide)
ch = 4;
tr = 1:10;
fr = 23; 
newFigure(1)
subplot 211
plot(squeeze(dataSlide(:,tr,ch,fr)))
subplot 212 
plot(squeeze(sin(dataSlide(:,tr,ch,fr))))


newFigure(2)
subplot 211
plot(squeeze(dataUnwrap(:,tr,ch,fr)))
subplot 212 
plot(squeeze(sin(dataUnwrap(:,tr,ch,fr))))


newFigure(3)
subplot 211
plot(squeeze(data0(:,tr,ch,fr)))
subplot 212 
plot(squeeze(sin(data0(:,tr,ch,fr))))

newFigure(4)
time = linspace(0,ep_size,size(datOut,1));
freq = foi;
dat  = squeeze(datOut(:,4,:));
imagesc(time,freq,dat')
set(gca,'clim',[0 .3])


%%
centClim     = @(x) set(x,'clim',max(abs(get(x,'clim')))*[-1 1]);

newFigure
ch   = 4;
max_val = .8;
dat1 = cat(4,DAT{1});
dat2 = cat(4,DAT{2});

dat1 = squeeze(mean(mean(dat1(:,ch,:,:),2),4));
dat2 = squeeze(mean(mean(dat2(:,ch,:,:),2),4));

% close all
time = linspace(0,ep_size,size(datOut,1));
freq = foi;
subplot 131
dat  = squeeze(dat1);
imagesc(time,freq,dat')
set(gca,'clim',[0 max_val])
title('ISOTONIC')
subplot 132
dat  = squeeze(dat2);
imagesc(time,freq,dat')
set(gca,'clim',[0 max_val])
title('DYNAMIC')
subplot 133
dat  = squeeze(dat1./dat2)-1;
imagesc(time,freq,dat')
centClim(gca)
title('CONTRAST')
