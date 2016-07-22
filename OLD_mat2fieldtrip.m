fs              = hdr.Fs;
cursor          = out.cursor;
target          = out.target;
samples         = out.samples(1):out.samples(end);
if ~isfield(out,'data')
    out.data           = double(ft_get_data(samples(1),samples(end)));
end

data.label             = hdr.label;
data.fsample           = fs;
data.time{tr}          = (samples-min(samples))/fs;
data.trial{tr}         = cat(1,out.data,cursor,target);
data.trialinfo(tr,1)   = k_trial;
data.sampleinfo(tr,:)  = out.samples;

if isfield(data,'trl')
     data.trl          = cat(1,data.trl,out.trl);
else data.trl          = out.trl; end






%{
% %%% create sampleinfo
% if isfield(data,'sampleinfo')
%     data.sampleinfo(end+1,1)  = samples(1);
%     data.sampleinfo(end+1,2)  = samples(1) + size(data.trial{1},2)-1;
% else
%     data.sampleinfo(:,1)     = samples(1);
%     data.sampleinfo(:,2)     = samples(1) + size(data.trial{1},2)-1;
% end
%}
