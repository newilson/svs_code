function f = NWman_phase(spec,ax_vals,fig_title)

spec_orig = spec;

if ndims(spec)>2
    errordlg('input cannot be more than 2D')
end

isvec = false;
if isvector(spec)==1
    isvec = true;
    spec = spec(:);
end
si = size(spec); % update

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
ax = axes('Parent',f,'position',[.1 .3 .7 .6]);
h = plot(ax,ax_vals,squeeze(real(spec(:,1)))); 
yl = [1.1*min(real(spec(:))), 0.9*max(real(spec(:)))];
yl = sort(yl); % in case min>max
ylim(ax,yl);
xlim(ax,[min(ax_vals) max(ax_vals)]);
set(ax,'XDir','reverse');


if ~isvector(spec)
    s1 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.9 .1 .05 .55],...
        'value',1,'min',1,'max',si(2),...
        'sliderstep',[1 1]/(si(2)-1),'callback',@nextslice);
    t1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.9 .05 .05 .05],...
        'string',num2str(1));
end
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

updateStrings(yl)

te1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .70 .05 .05],...
    'string','Min','HorizontalAlignment','right','fontweight','bold');

te2 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .65 .05 .05],...
    'string','Max','HorizontalAlignment','right','fontweight','bold');

s2 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.1 .15 .55 .05],...
        'value',0,'min',-179.9,'max',180.0,...
        'sliderstep',[1/360 10/360],'callback',@updatePC0sl);
    
s3 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.1 .05 .55 .05],...
        'value',0,'min',-1079.9,'max',1080.0,...
        'sliderstep',[1/500 5/100],'callback',@updatePC1sl);
    
te3 = uicontrol('Parent',f,'style','text','units','normalized','position',[.05 .15 .05 .05],...
    'string','PC0','HorizontalAlignment','left','fontweight','bold');

te4 = uicontrol('Parent',f,'style','text','units','normalized','position',[.05 .05 .05 .05],...
    'string','PC1','HorizontalAlignment','left','fontweight','bold');

e3 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.7 .15 .1 .05],...
    'callback',@updatePC0ed);
set(e3,'string',num2str(0));

e4 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.7 .05 .1 .05],...
    'callback',@updatePC1ed);
set(e4,'string',num2str(0));

te5 = uicontrol('Parent',f,'style','text','units','normalized','position',[.1 .93 .05 .05],...
    'string','Pivot','HorizontalAlignment','left','fontweight','bold');

e5 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.15 .93 .1 .05],...
    'callback',@updatePCpiv);
set(e5,'string',num2str(ax_vals(end/2+1)));

b3 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.5 .93 .1 .05],...
    'callback',@save_var,'string','Save','fontweight','bold');

if isvec
    slice = 1;
end

set(f,'visible','on','toolbar','figure')

    function nextslice(source,callbackdata)
        slice = round(get(source,'value'));
        set(h,'ydata',squeeze(real(spec(:,slice))))
        set(t1,'string',num2str(slice))
        set(s1,'value',slice)
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

    function updatePC0sl(source,callbackdata) % updated slider
        set(e3,'string',num2str(get(s2,'value')));
        doPC;
    end

    function updatePC1sl(source,callbackdata) % updated slider
        set(e4,'string',num2str(get(s3,'value')));
        doPC;
    end

    function updatePC0ed(source,callbackdata) % updated edit box
        set(s2,'value',str2double(get(e3,'string')));
        doPC;
    end

    function updatePC1ed(source,callbackdata) % updated edit box
        set(s3,'value',str2double(get(e4,'string')));
        doPC;
    end

    function updatePCpiv(source,callbackdata) % updated pivot point
        doPC;
    end

    function doPC
        pivot = str2double(get(e5,'string'));
        pc0 = str2double(get(e3,'string'));
        pc1 = str2double(get(e4,'string'));
        linph = 2*pc1/abs(ax_vals(end)-ax_vals(1)) * (ax_vals-pivot);
        phasevec = exp(-1i*pi/180*(pc0 + linph));
        phasearr = repmat(phasevec,[1,si(2)]);
        spec = spec_orig .* phasearr;
        if ~isvec
            slice = get(s1,'value');
        end
        set(h,'ydata',squeeze(real(spec(:,slice))))
    end

    function save_var(source,callbackdata)
        assignin('caller','out',spec);
        parsPC.pivot = str2double(get(e5,'string'));
        parsPC.pc0 = str2double(get(e3,'string'));
        parsPC.pc1 = str2double(get(e4,'string')); 
        assignin('caller','parsPC',parsPC);
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