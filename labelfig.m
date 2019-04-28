function labelfig(fighand,figname,figdir,unlabel)
% labels figure with filename

if nargin < 3
    figdir=pwd;
end
if nargin < 4
    unlabel=0; 
end
if strcmp(figdir,['.' filesep])
    figdir=pwd;  % replace by full path name
end

figure(fighand)

if ~unlabel
    if isempty(figdir)
        flab=figname;
    else
        if strcmp(figdir(end),filesep)
            flab=[figdir figname];
        else
            flab=[figdir filesep figname];
        end
        maxfolders=4;
        sepix=regexp(flab,filesep);
        if ~isempty(sepix)
            totsep=length(sepix);
            ixfirst=sepix(totsep-min(maxfolders,totsep-1));
            flab=['..' flab(ixfirst:end)];
        end
    end
    
    
    subplot('position',[0 0 1 0.05])
    axis([0 1 0 1])
    delete(findobj(fighand,'Tag','directoryLabel'))
    hdir=text(0,0.05,flab,'verticalalignment','bottom','fontsize',6,...
        'interpreter','none');
    set(hdir,'Tag','directoryLabel')
    delete(findobj(fighand,'Tag','dateLabel'))
    hdate=text(1,0.05,date,'horizontalalignment','right',...
        'verticalalignment','bottom','fontsize',6);
    set(hdate,'Tag','dateLabel')
    axis off
else  % remove existing labels
    delete(findobj(fighand,'Tag','directoryLabel'))
    delete(findobj(fighand,'Tag','dateLabel'))
end
