function f = NWsvdTS(spec,ax_vals,fig_title)

spec = double(spec); % must be double

spec_orig = spec;

if ndims(spec)>2
    errordlg('input cannot be more than 2D')
end

if isvector(spec)==1
    errordlg('input must be time series of spectra')
end

si = size(spec);
if si(1)<si(2)
    warndlg(['spectral points = ' num2str(si(1)) ', and time points = ' num2str(si(2))]);
end

if nargin<2 || isempty(ax_vals) || ~isequal(length(ax_vals),si(1))
    ax_vals = 1:si(1);
end
ax_vals = ax_vals(:);

f = figure('position',[200 50 800 580]);
if nargin==3 && ~isempty(fig_title)
    if iscell(fig_title)
        set(f,'Name',fig_title{1},'NumberTitle','off')
    else
        set(f,'Name',fig_title,'NumberTitle','off');
    end
else
    fig_title = [];
end

s1 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.9 .1 .05 .55],...
    'value',1,'min',1,'max',si(2),...
    'sliderstep',[1 1]/(si(2)-1),'callback',@nextslice);
t1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.9 .05 .05 .05],...
    'string',num2str(1));
    
b0 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .93 .07 .05],...
    'callback',@printax,'string','Print','fontweight','bold');

b1 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .85 .05 .05],...
    'callback',@ylim_2,'string','/2');

b2 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.9 .80 .05 .05],...
    'callback',@ylimX2,'string','x2');

e1 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.9 .70 .07 .05],...
    'callback',@ylimMin);

e2 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.9 .65 .07 .05],...
    'callback',@ylimMax);


te1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .70 .05 .05],...
    'string','Min','HorizontalAlignment','right','fontweight','bold');

te2 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .65 .05 .05],...
    'string','Max','HorizontalAlignment','right','fontweight','bold');  
    

te5 = uicontrol('Parent',f,'style','text','units','normalized','position',[.1 .93 .1 .05],...
    'string','Ready...','HorizontalAlignment','left','fontweight','bold');


b3 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.5 .93 .1 .05],...
    'callback',@save_var,'string','Save','fontweight','bold');

te7 = uicontrol('Parent',f,'style','text','units','normalized','position',[.27 .15 .1 .05],...
    'string','Singular Values','HorizontalAlignment','left','fontweight','bold');

def_sv = 5;
e4 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.35 .16 .05 .04],'callback',@update);
set(e4,'string',num2str(def_sv));

set(f,'visible','on','toolbar','figure')

bg = uibuttongroup('Visible','on','position',[.65 .92 .15 .1],'BorderType','none','SelectionChangedFcn',@updateRadio);
r1 = uicontrol('Parent',bg,'style','radiobutton','units','normalized','position',[.05 .2 .5 .5],...
    'string','Original');
r2 = uicontrol('Parent',bg,'style','radiobutton','units','normalized','position',[.6 .2 .4 .5],...
    'string','Diff');
r1.Value = true; % default is Original

% initial denoising
[U,S,V] = svds(spec_orig,def_sv,'largest','MaxIterations',250);
spec = U*S*V';

ax = axes('Parent',f,'position',[.1 .3 .7 .6]);
h = plot(ax,ax_vals,squeeze(real(spec_orig(:,1))),'k-',ax_vals,squeeze(real(spec(:,1))),'r-'); 
set(h(1),'linewidth',0.5,'color',[0.5 0.5 0.5]); % gray
set(h(2),'linewidth',1.0);
yl = [1.1*min(real(spec_orig(:))), 0.9*max(real(spec_orig(:)))];
yl = sort(yl); % in case min>max
ylim(ax,yl);
xlim(ax,[min(ax_vals) max(ax_vals)]);
set(ax,'XDir','reverse');


updateStrings(yl)

    function nextslice(source,callbackdata)
        slice = round(get(source,'value'));
        set(t1,'string',num2str(slice))
        set(s1,'value',slice)
        doPlot;
        if iscell(fig_title)
            set(f,'Name',fig_title{slice})
        end
    end

    function ylim_2(source,callbackdata)
        yl = ylim(ax)/2;
        ylim(ax,yl);
        updateStrings(yl)
    end
    
    function ylimX2(source,callbackdata)
        yl = ylim(ax)*2;
        ylim(ax,yl);
        updateStrings(yl)
    end

    function ylimMax(source,callbackdata)
        yl = ylim(ax);
        if str2double(get(e2,'string'))>yl(1)
            yl(2) = str2double(get(e2,'string'));
            ylim(ax,yl);
        else
            updateStrings(yl)
        end
    end

    function ylimMin(source,callbackdata)
        yl = ylim(ax);
        if str2double(get(e1,'string'))<yl(2)
            yl(1) = str2double(get(e1,'string'));
            ylim(ax,yl);
        else
            updateStrings(yl)
        end
    end

    function updateStrings(yl)
        set(e1,'string',num2str(yl(1),3))
        set(e2,'string',num2str(yl(2),3))
    end

    function update(source,callbackdata)
        set(te5,'string','filtering...');
        nsv = round(str2double(e4.String));
        if nsv<1
            nsv = 1;
        elseif nsv>si(1)
            nsv = si(1);
        end
        e4.String = num2str(nsv);
        
        [funU,funS,funV] = svds(spec_orig,nsv,'largest','MaxIterations',250);
        spec = funU*funS*funV';

        doPlot;
        
        set(te5,'string','DONE');
    end  

    function updateRadio(source,callbackdata)
        doPlot;
    end

    function doPlot
        slice = get(s1,'value');
        if r1.Value
            set(h(1),'ydata',squeeze(real(spec_orig(:,slice))))
        elseif r2.Value
            set(h(1),'ydata',squeeze(real(spec_orig(:,slice)-spec(:,slice))))
        end
        set(h(2),'ydata',squeeze(real(spec(:,slice))))
    end

    function save_var(source,callbackdata)
        assignin('caller','out',spec);
    end

    function printax(source,callbackdata)
        prompt = {'FullName (no extension)','Format (vector: eps,pdf / bitmap: tiff,png,bmp,jpeg)','DPI','Renderer (painters or opengl)'};
        defvals = {fullfile(pwd,'myFigure'),'eps','300','painters'};
        nlines = 1;
        vals = inputdlg(prompt,'Options',nlines,defvals);
        if ~isempty(vals)
            fname = vals{1};
            if strcmp(vals{2},'eps') && ~strcmp(cmap,'bone')
                form = '-depsc';
            else
                form = ['-d' vals{2}];
            end
            res = ['-r' vals{3}];
            rend = ['-' vals{4}]; 
            print(f,fname,form,res,rend,'-noui')
        end
    end
  
end

