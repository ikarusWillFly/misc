function plot_coh_arrows(cfg,statSensor)
hold on
%%% input checking
field = 'alpha';        value = 0.025;
if ~isfield(cfg,field), cfg.(field) = value; end

lay    = cfg.hdr.layout;
cmap   = colormap;
cmap1  = cmap(1,:);
cmap2  = cmap(end,:);
sig    = find(statSensor.prob <= cfg.alpha);

for j = 1 : numel(sig)
    ind = sig(j);
    chans  = statSensor.labelcmb(ind,:);
    chInd1 = lay.pos(strcmpi(lay.label,chans(1)),:);
    chInd2 = lay.pos(strcmpi(lay.label,chans(2)),:);
    x = [chInd1(1),chInd2(1)];
    y = [chInd1(2),chInd2(2)];
    delta = chInd2-chInd1;
    
    tmp = statSensor.stat(ind);
    
    if statSensor.stat(ind) > 0
        color = cmap1;
    else
        color = cmap2;
    end
    patchline(x,y,'linewidth',abs(tmp),'edgecolor',color,'edgealpha',0.3)
%     plot(x,y,'color',color,'linewidth',abs(tmp),'markeredgecolor','k')
%     plot(x(end),y(end),'d','color',color,'markerfacecolor',color)
end

end
