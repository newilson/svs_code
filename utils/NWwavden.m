function f = NWwavden(spec,ax_vals,fig_title)

% check for files
if ~exist('wmaxlev')
    warndlg('cannot find wmaxlev');
    return;
end
newver = true;
if ~exist('wdenoise')
    if ~exist('wden')
        warndlg('cannot find wdenoise or wden');
        return;
    end
    newver = false;
end

spec = double(spec); % must be double for wdenoise

spec_orig = spec;

if ndims(spec)>2
    errordlg('input cannot be more than 2D')
end

isRe = isreal(spec);

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

lev_max = max(floor(log2(si(1))),wmaxlev(si(1),'sym4'));
lev_def = min(floor(log2(si(1))),wmaxlev(si(1),'sym4'));

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


te1 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .70 .05 .05],...
    'string','Min','HorizontalAlignment','right','fontweight','bold');

te2 = uicontrol('Parent',f,'style','text','units','normalized','position',[.84 .65 .05 .05],...
    'string','Max','HorizontalAlignment','right','fontweight','bold');


% s2 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.1 .15 .55 .05],...
%         'value',lev_def,'min',1,'max',lev_max,...
%         'sliderstep',[1 1]/(lev_max-1),'callback',@updateLev_sl);    
    
te3 = uicontrol('Parent',f,'style','text','units','normalized','position',[.5 .15 .05 .05],...
    'string','Level','HorizontalAlignment','left','fontweight','bold');

e3 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.55 .16 .05 .04],...
    'callback',@updateLev_ed);
set(e3,'string',num2str(lev_def));


te5 = uicontrol('Parent',f,'style','text','units','normalized','position',[.1 .93 .1 .05],...
    'string','Ready...','HorizontalAlignment','left','fontweight','bold');


b3 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.5 .93 .1 .05],...
    'callback',@save_var,'string','Save','fontweight','bold');

te6 = uicontrol('Parent',f,'style','text','units','normalized','position',[.1 .15 .1 .05],...
    'string','Method','HorizontalAlignment','left','fontweight','bold');

pop1 = uicontrol('Parent',f,'style','popupmenu','units','normalized','position',[.17 .1 .1 .1],'callback',@updateMethod);
if newver
    pop1.String = {'Bayes','Minimax','SURE','UniversalThreshold'};
    pop1.Value = 1;
else
    pop1.String = {'sqtwolog','minimaxi','heursure'};
    pop1.Value = 1;
end

te7 = uicontrol('Parent',f,'style','text','units','normalized','position',[.3 .15 .1 .05],...
    'string','Rule','HorizontalAlignment','left','fontweight','bold');

pop2 = uicontrol('Parent',f,'style','popupmenu','units','normalized','position',[.35 .1 .1 .1],'callback',@updateRule);
if newver
    pop2.String = {'Hard','Soft','Median','Mean'};
    pop2.Value = 3;
else
    pop2.String = {'h','s'};
    pop2.Value = 1;
end

if isvec
    slice = 1;
end

set(f,'visible','on','toolbar','figure')

bg = uibuttongroup('Visible','on','position',[.65 .92 .15 .1],'BorderType','none','SelectionChangedFcn',@updateRadio);
r1 = uicontrol('Parent',bg,'style','radiobutton','units','normalized','position',[.05 .2 .5 .5],...
    'string','Original');
r2 = uicontrol('Parent',bg,'style','radiobutton','units','normalized','position',[.6 .2 .4 .5],...
    'string','Diff');
r1.Value = true; % default is Original


% initial denoising
for ii=1:si(2)
    if newver
        tempR = wdenoise(real(spec_orig(:,ii)),lev_def,'Wavelet','sym4','DenoisingMethod','Bayes','ThresholdRule','Median','NoiseEstimate','LevelIndependent');
        if ~isRe
            tempI = wdenoise(imag(spec_orig(:,ii)),lev_def,'Wavelet','sym4','DenoisingMethod','Bayes','ThresholdRule','Median','NoiseEstimate','LevelIndependent');
            spec(:,ii) = complex(tempR,tempI);
        else
            spec(:,ii) = tempR;
        end
    else
        tempR = wden(real(spec_orig(:,ii)),'sqtwolog','h','sln',lev_def,'sym4');
        if ~isRe
            tempI = wden(imag(spec_orig(:,ii)),'sqtwolog','h','sln',lev_def,'sym4');
            spec(:,ii) = complex(tempR,tempI);
        else
            spec(:,ii) = tempR;
        end
    end
end
% h = plot(ax,ax_vals,squeeze(real(spec(:,1))),'r-');

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

%     function updateLev_sl(source,callbackdata) % updated slider
%         set(te5,'string','filtering...');
%         set(e3,'string',num2str(get(s2,'value')));
%         doWavDen;
%         set(te5,'string','DONE');
%     end

    function updateLev_ed(source,callbackdata) % updated edit box
        cur_lev = round(str2double(get(e3,'string')));
        if cur_lev < 1
            cur_lev = 1;
        elseif cur_lev > lev_max
            cur_lev = lev_max;
        end
        set(e3,'string',num2str(cur_lev)); % round
%         set(s2,'value',str2double(get(e3,'string')));
        doWavDen;
    end

    function updateMethod(source,callbackdata)
        cur_ruleval = pop2.Value;
        if cur_ruleval>2
            pop2.Value = 1;
        end
        doWavDen;
    end

    function updateRule(source,callbackdata)
        doWavDen;
    end

    function doWavDen
        set(te5,'string','filtering...');
        cur_lev = str2double(get(e3,'string'));
        cur_methval = pop1.Value;
        strmeth = pop1.String;
        cur_meth = strmeth{cur_methval};
        cur_ruleval = pop2.Value;
        strrule = pop2.String;
        cur_rule = strrule{cur_ruleval};
        for funind=1:si(2)
            if newver
                funtempR = wdenoise(real(spec_orig(:,funind)),cur_lev,'Wavelet','sym4','DenoisingMethod',cur_meth,'ThresholdRule',cur_rule,'NoiseEstimate','LevelIndependent');
                if ~isRe
                    funtempI = wdenoise(imag(spec_orig(:,funind)),cur_lev,'Wavelet','sym4','DenoisingMethod',cur_meth,'ThresholdRule',cur_rule,'NoiseEstimate','LevelIndependent');
                    spec(:,funind) = complex(funtempR,funtempI);
                else
                    spec(:,funind) = funtempR;
                end
            else
                funtempR = wden(real(spec_orig(:,ii)),cur_meth,cur_rule,'sln',cur_lev,'sym4');
                if ~isRe
                    funtempI = wden(imag(spec_orig(:,ii)),cur_meth,cur_rule,'sln',cur_lev,'sym4');
                    spec(:,ii) = complex(funtempR,funtempI);
                else
                    spec(:,ii) = funtempR;
                end
            end
        end
        doPlot;
        set(te5,'string','DONE');
    end

    function updateRadio(source,callbackdata)
        doPlot;
    end

    function doPlot
        if ~isvec
            slice = get(s1,'value');
        end
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