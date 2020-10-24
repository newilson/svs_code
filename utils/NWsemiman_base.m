function f = NWsemiman_base(spec,ax_vals,fig_title)

spec = double(spec);
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

s2 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.1 .15 .45 .05],...
        'value',-2,'min',-5,'max',3,...
        'sliderstep',[0.1 0.5]/8,'callback',@updatePsl);
    
s3 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.1 .05 .45 .05],...
        'value',5,'min',-1,'max',12,...
        'sliderstep',[0.1 0.5]/13,'callback',@updateLsl);
    
te3 = uicontrol('Parent',f,'style','text','units','normalized','position',[.05 .15 .05 .05],...
    'string','log p','HorizontalAlignment','left','fontweight','bold');

te4 = uicontrol('Parent',f,'style','text','units','normalized','position',[.05 .05 .05 .05],...
    'string','log L','HorizontalAlignment','left','fontweight','bold');

e3 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.6 .15 .05 .05],...
    'callback',@updatePed);
set(e3,'string',num2str(-2));

e4 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.6 .05 .05 .05],...
    'callback',@updateLed);
set(e4,'string',num2str(5));

% te5 = uicontrol('Parent',f,'style','text','units','normalized','position',[.1 .93 .05 .05],...
%     'string','Pivot','HorizontalAlignment','left','fontweight','bold');
% 
% e5 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.15 .93 .1 .05],...
%     'callback',@updatePCpiv);
% set(e5,'string',num2str(ax_vals(end/2+1)));

b3 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.5 .93 .1 .05],...
    'callback',@save_var,'string','Save','fontweight','bold');

minsize = max(diff(ax_vals(:)));
maxsize = (max(ax_vals(:))-min(ax_vals(:)))/2 - minsize;
defsize = min(maxsize/5,200);
s4 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.12 .15 .45 .05],...
        'value',defsize,'min',minsize,'max',maxsize,...
        'sliderstep',[0.02 0.10],'callback',@updateStepsl);
    
s5 = uicontrol('Parent',f,'Style','slider','units','normalized','Position',[.12 .05 .45 .05],...
        'value',defsize,'min',minsize,'max',maxsize,...
        'sliderstep',[0.02 0.10],'callback',@updateWinsl);
    
te5 = uicontrol('Parent',f,'style','text','units','normalized','position',[.05 .15 .055 .05],...
    'string','Step Size','HorizontalAlignment','left','fontweight','bold');

te6 = uicontrol('Parent',f,'style','text','units','normalized','position',[.05 .05 .055 .05],...
    'string','Window Size','HorizontalAlignment','left','fontweight','bold');

e5 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.6 .15 .05 .05],...
    'callback',@updateSteped);
set(e5,'string',num2str(s4.Value));

e6 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.6 .05 .05 .05],...
    'callback',@updateWined);
set(e6,'string',num2str(s5.Value));

te7 = uicontrol('Parent',f,'style','text','units','normalized','position',[.7 .15 .1 .05],...
    'string','Regression Method','HorizontalAlignment','left','fontweight','bold');

pop1 = uicontrol('Parent',f,'style','popupmenu','units','normalized','position',[.8 .15 .07 .05],...
    'string',{'pchip','linear','spline'},'value',1,'callback',@updateMethod);

te8 = uicontrol('Parent',f,'style','text','units','normalized','position',[.7 .08 .1 .05],...
    'string','Estimation Method','HorizontalAlignment','left','fontweight','bold');

pop2 = uicontrol('Parent',f,'style','popupmenu','units','normalized','position',[.8 .08 .07 .05],...
    'string',{'quantile','em'},'value',1,'callback',@updateMethod);

te9 = uicontrol('Parent',f,'style','text','units','normalized','position',[.7 .02 .1 .05],...
    'string','Quantile','HorizontalAlignment','left','fontweight','bold');

e7 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.8 .04 .04 .03],...
    'callback',@updateMethodQu);
set(e7,'string',0.10);

bg = uibuttongroup('Visible','on','position',[.05 .9 .4 .1],'BorderType','none','SelectionChangedFcn',@updateMethod);
r1 = uicontrol('Parent',bg,'style','radiobutton','units','normalized','position',[.05 .2 .5 .5],...
    'string','Asymmetric Least Squares');
r2 = uicontrol('Parent',bg,'style','radiobutton','units','normalized','position',[.65 .2 .4 .5],...
    'string','Back Adjust');
r1.Value = true; % default is ALS

updateVisibility;
doBaselineALS;

if ~exist('msbackadj') % check for file
    bg.Visible = 'off';
end

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

    function updatePsl(source,callbackdata) % updated slider
        set(e3,'string',num2str(get(s2,'value')));
        doBaselineALS;
    end

    function updateLsl(source,callbackdata) % updated slider
        set(e4,'string',num2str(get(s3,'value')));
        doBaselineALS;
    end

    function updatePed(source,callbackdata) % updated edit box
        set(s2,'value',str2double(get(e3,'string')));
        doBaselineALS;
    end

    function updateLed(source,callbackdata) % updated edit box
        set(s3,'value',str2double(get(e4,'string')));
        doBaselineALS;
    end

    function updateStepsl(source,callbackdata) % updated slider
        set(e5,'string',num2str(get(s4,'value')));
        doBaselineBA;
    end

    function updateWinsl(source,callbackdata) % updated slider
        set(e6,'string',num2str(get(s5,'value')));
        doBaselineBA;
    end

    function updateSteped(source,callbackdata) % updated edit box
        set(s4,'value',str2double(get(e5,'string')));
        doBaselineBA;
    end

    function updateWined(source,callbackdata) % updated edit box
        set(s5,'value',str2double(get(e6,'string')));
        doBaselineBA;
    end

    function updateMethodQu(source,callbackdata)
        val = str2double(e7.String);
        if val<0.01
            e7.String = num2str(0.01);
        elseif val>1
            e7.String = num2str(1);
        end
        doBaselineBA; % only this option
    end

    function updateMethod(source,callbackdata) % updated radio button
        updateVisibility;
        if r1.Value
            doBaselineALS;
        elseif r2.Value
            doBaselineBA;
        end
    end

    function updateVisibility
        if r1.Value
            
            % ALS options
            te3.Visible = 'on';
            te4.Visible = 'on';
            e3.Visible = 'on';
            e4.Visible = 'on';
            s2.Visible = 'on';
            s3.Visible = 'on';
            
            % BA options
            te5.Visible = 'off';
            te6.Visible = 'off';
            te7.Visible = 'off';
            te8.Visible = 'off';
            te9.Visible = 'off';
            e5.Visible = 'off';
            e6.Visible = 'off';
            e7.Visible = 'off';
            s4.Visible = 'off';
            s5.Visible = 'off';
            pop1.Visible = 'off';
            pop2.Visible = 'off';
            
        elseif r2.Value
            
            % ALS options
            te3.Visible = 'off';
            te4.Visible = 'off';
            e3.Visible = 'off';
            e4.Visible = 'off';
            s2.Visible = 'off';
            s3.Visible = 'off';
            
            % BA options
            te5.Visible = 'on';
            te6.Visible = 'on';
            te7.Visible = 'on';
            te8.Visible = 'on';
            te9.Visible = 'on';
            e5.Visible = 'on';
            e6.Visible = 'on';
            e7.Visible = 'on';
            s4.Visible = 'on';
            s5.Visible = 'on';
            pop1.Visible = 'on';
            pop2.Visible = 'on';
        end
    end

    function doBaselineALS
        p = 10^str2double(get(e3,'string'));
        lambda = 10^str2double(get(e4,'string'));
        spec = 0*spec_orig;
        for ii=1:si(2)
            baseline = baselinecalc(real(spec_orig(:,ii)),lambda,p);
            spec(:,ii) = complex(real(spec_orig(:,ii)) - baseline, imag(spec_orig(:,ii)) - baseline); % same baseline subtraction for real/imag
        end
        if ~isvec
            slice = get(s1,'value');
        else
            slice = 1;
        end
        set(h,'ydata',squeeze(real(spec(:,slice))))
    end

    function doBaselineBA
        regval = pop1.Value;
        regstr = pop1.String;
        regmeth = regstr{regval};
        estval = pop2.Value;
        eststr = pop2.String;
        estmeth = eststr{estval};
        stepsize = s4.Value;
        winsize = s5.Value;
        quantile = str2double(e7.String);
        spec = 0*spec_orig;
        for ii=1:si(2)
            spec(:,ii) = msbackadj(ax_vals(:),real(spec_orig(:,ii)),'RegressionMethod',regmeth,'EstimationMethod',estmeth,...
                'StepSize',stepsize,'WindowSize',winsize,'QuantileValue',quantile,'SmoothMethod','none','ShowPlot',false);
            if ~isRe
                baseline = real(spec_orig(:,ii))-spec(:,ii);
                spec(:,ii) = complex(spec(:,ii),imag(spec_orig(:,ii))-baseline); % same baseline subtraction for real/imag
            end
        end
        if ~isvec
            slice = get(s1,'value');
        end
        set(h,'ydata',squeeze(real(spec(:,slice))));
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

function z = baselinecalc(y, lambda, p)
% Estimate baseline with asymmetric least squares

    m = length(y);
    D = diff(speye(m), 2);
    w = ones(m, 1);
    
    for it = 1:6
        W = spdiags(w, 0, m, m);
        C = chol(W + lambda * (D' * D));
        z = C \ (C' \ (w .* y));
        w = p * (y > z) + (1 - p) * (y < z);
    end

end