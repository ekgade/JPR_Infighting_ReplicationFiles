function net2dplot(x,y,adjmat,nodcol,groups,titlestr,axlabs,thresh)
% net2dplot(x,y,adjmat,nodcol,groups,titlestr,axlabs,thresh)
% Plots nodes on 2-d plane with lines showing network connections
% Inputs:
% x,y: x and y coordinates of nodes
% adjmat: Network adjacency matrix (symmetrized in function
% if not already so). Set to empty to not plot links.
% nodcol: array containing color triples for nodes
% groups: names corresponding to nodes
% titlestr: title for plot (default is none)
% axlabs: 2-d cell containing labels for x and y axes (default 
% is none)
% thresh: tie strength threshold for plotting links

if nargin<6
    titlestr='';
end
if nargin<7 || isempty(axlabs)
    axlabs={'',''};
end
if nargin<8
    thresh=0;
end


N=length(x);


figure
hold on
if ~isempty(adjmat)
    im_mut=(adjmat+adjmat')/2;  % symmetrize
    lnsty='-';
    lw=im_mut/mean(mean(im_mut(im_mut>0))); % normalize by mean link strength
    for m=1:N
        for n=1:m-1
            if im_mut(m,n)>thresh
                line([x(m) x(n)],[y(m) y(n)],...
                    'color','k','linestyle',lnsty,'linewidth',lw(m,n));
            end
        end
    end
end

for k=1:N
    line(x(k),y(k),'marker','o','markerfacecolor',nodcol(k,:),...
        'markeredgecolor','k','markersize',10)   
end
ax=axis;
for k=1:N
    text(x(k)+0.02*(ax(2)-ax(1)),y(k),groups{k}...
        ,'fontweight','bold','fontsize',10)
end

hold off
set(gca,'fontsize',12)
xlabel(axlabs{1});ylabel(axlabs{2})
title(titlestr)

