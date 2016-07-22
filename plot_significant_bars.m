function h = plot_significant_bars(significance_matrix,color)
if nargin < 2
    color = 'k';
end
ca = gca;
ylims  = max(ca.YLim);
tmp    = sign(ylims);
xdata  = ca.XTick;
if any(significance_matrix(:))
    [sig(:,1),sig(:,2)] = find(significance_matrix);
    sig                 = unique(sort(sig,2),'rows');
    for k    = 1 : size(sig,1)
        tmp2 = xdata(sig(k,:));
        tmp2 = tmp2+[+1 -1]*.03;
        x    = ([tmp2(1),tmp2,tmp2(end)]);
        %         y    = + zeros(1,numel(x)) + ca.YLim(2) * [.99 1 1 .99] + ca.YLim(2)*k/40;
        y    = ca.YLim(2) * (ones(1,4)+ range(ca.YLim)*[.05 0 0 .05].*-tmp);% - ca.YLim(2)*k/40.*-tmp;
        h(k) = plot(x,y,'-','color',color,'linewidth',2);
        %     text(mean(x),max(y)*1.0,'*','fontSize',30)
        cumY(k,:) = y;
    end
    if  tmp>0
     ylim([min(ylim) max(cumY(:))+[.01*tmp]])

    else
    ylim([min(ylim) max(cumY(:))+[.01*-tmp]])
    end
    
end
end