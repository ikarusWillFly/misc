
function fit = IO_curve(data,x,intX)
rep      = size(data,1);
sigm     = @(param,xval) (param(1)./(1+10.^((param(2)-xval)*param(3)))) +param(4);

meanD    = nanmean(data);

params                          = [max(meanD),x(find(diff(meanD) == max(diff(meanD)))+1),.5,min(meanD)];
[fit.BETA,fit.R,~,fit.COVB,~]   = nlinfit(reshape(repmat(x,rep,1),[],1),data(:),sigm,params);
[fit.yPRED, ~]                  = nlpredci(sigm,intX,fit.BETA,fit.R,'covar',fit.COVB);
end