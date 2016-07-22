% clear all
% datIn = randn(256*2,10,50,20)*pi/2;

datIn = unwrap(angle(datIn));

%%%  CUT TRIALS IN SLIDING WINDOWS
cfg = [];
cfg.time_dimension  = 1;              % time dimension
cfg.trial_dimension = 4;              % trial dimension
cfg.epoch_size      = round(.100*256);% desired epoch sizes (in samples)
cfg.sliding         = 25;             % desired sliding size (in samples)
cfg.keeptrials      = false;          % flag to keep the trials separated
cfg.feedback        = true;           % flag to receive feedback with datadimensions
datOut              = slidingEpochs(cfg,datIn); clear datIn

%%% MANIPULATE THE DATA
datOut              = bsxfun(@minus,datOut,datOut(1,:,:,:)); % baseline the epochs to the first phase
datOut              = permute(datOut,[2 1 3 4]);             % get the repetitions dimension first 
datOut              = squeeze(std(datOut));                  % perform the standard deviation 
datOut              = (2*datOut/(2*pi))*100;                 % calculate the percentage of a cycle equal to 2 std of the phases
plot(datOut(:,:,1))                                          % plot it