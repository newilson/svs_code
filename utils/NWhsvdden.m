function f = NWhsvdden(fid,ax_vals,bw,fig_title)

fid = double(fid); % must be double

fid_orig = fid;

if ndims(fid)>2
    errordlg('input cannot be more than 2D')
end

isvec = false;
if isvector(fid)==1
    isvec = true;
    fid = fid(:);
end
si = size(fid); % update

if nargin<2 || isempty(ax_vals) || ~isequal(length(ax_vals),si(1))
    ax_vals = 1:si(1);
end
ax_vals = ax_vals(:);

if nargin<3 || isempty(bw) || numel(bw)>1, bw = 1; end
time = (0:si(1)-1)*1/bw;

f = figure('position',[200 50 800 580]);
if nargin==4 && ~isempty(fig_title)
    if iscell(fig_title)
        set(f,'Name',fig_title{1},'NumberTitle','off')
    else
        set(f,'Name',fig_title,'NumberTitle','off');
    end
else
    fig_title = [];
end

if ~isvector(fid)
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
    
te3 = uicontrol('Parent',f,'style','text','units','normalized','position',[.5 .15 .05 .05],...
    'string','LB','HorizontalAlignment','left','fontweight','bold');

e3 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.55 .16 .05 .04],...
    'callback',@update);
set(e3,'string',num2str(0)); % default - no line broadening


te5 = uicontrol('Parent',f,'style','text','units','normalized','position',[.1 .93 .1 .05],...
    'string','Ready...','HorizontalAlignment','left','fontweight','bold');


b3 = uicontrol('Parent',f,'style','pushbutton','units','normalized','position',[.5 .93 .1 .05],...
    'callback',@save_var,'string','Save','fontweight','bold');

te7 = uicontrol('Parent',f,'style','text','units','normalized','position',[.27 .15 .1 .05],...
    'string','Singular Values','HorizontalAlignment','left','fontweight','bold');

def_sv = 5;
e4 = uicontrol('Parent',f,'style','edit','units','normalized','position',[.35 .16 .05 .04],'callback',@update);
set(e4,'string',num2str(def_sv));

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

% initial denoising - no line broadening
ncad = 2;
for ii=1:si(2)
    fid(:,ii) = NWhsvd(fid_orig(:,ii),def_sv,ncad);
end

spec_orig = fftshift(fft(fid_orig,[],1),1);
spec = fftshift(fft(fid,[],1),1);
ax = axes('Parent',f,'position',[.1 .3 .7 .6]);
h = plot(ax,ax_vals,squeeze(real(spec_orig(:,1))),'k-',ax_vals,squeeze(real(spec(:,1))),'r-'); 
set(h(1),'linewidth',0.5,'color',[0.5 0.5 0.5]); % gray
set(h(2),'linewidth',1.0);
yl = [1.1*min(real(spec_orig(:))), 0.9*max(real(spec_orig(:)))];
yl = sort(yl); % in case min>max
ylim(ax,yl);
xlim(ax,[min(ax_vals) max(ax_vals)]);

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
        
        lb = str2double(e3.String);
        fid = expFilter(time,lb,fid_orig);
        
        for funind=1:si(2)
            fid(:,funind) = NWhsvd(fid(:,funind),nsv,ncad);
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
        spec = fftshift(fft(fid,[],1),1);
        if r1.Value
            set(h(1),'ydata',squeeze(real(spec_orig(:,slice))))
        elseif r2.Value
            set(h(1),'ydata',squeeze(real(spec_orig(:,slice)-spec(:,slice))))
        end
        set(h(2),'ydata',squeeze(real(spec(:,slice))))
    end

    function save_var(source,callbackdata)
        assignin('caller','out',fid);
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

function [out,filt] = expFilter(t,lb,in)
%
% out = expFilter(t,lb,out)
%
% t is [npts x 1]
% lb is line broadening in Hz

if isequal(lb,0)
    filt = 1;
    out = in;
    return
end

if nargin>2
    npts = size(in,1);
    if ~isequal(length(t),npts)
        warning('number of points do not match')
        in = [];
    end
else
    in = [];
end

out = [];
filt = exp(-2*pi*t(:)*lb);

if ~isempty(in)
    si = size(in);
    in = reshape(in,npts,[]);
    filt2d = repmat(filt,[1 size(in,2)]);
    out = in .* filt2d;
    out = reshape(out,si);
end
end

function fid = NWhsvd(fid,nsv,ncad,ncol)
%
%

if ~isvector(fid), warning('not a vector'), return; end
doaug = false;
if mod(length(fid),2)==0, fid(end+1) = 0; doaug = true; end % make fid odd
npts = length(fid);

if nargin<2 || isempty(nsv), nsv = 10; end % number of singular values
if nargin<3 || isempty(ncad), ncad = 3; end % number of Cadzow iterations
if nargin<4 || isempty(ncol), ncol = ceil(npts/2); end % square Hankel matrix

nrow = npts - ncol + 1;

% construct Hankel matrix
H = zeros(nrow,ncol);
for ii=1:nrow
    H(ii,:) = fid(ii:ii+ncol-1);
end

% Cadzow iteration
for zz=1:ncad
    % truncated SVD
    [U,S,V] = svds(H,nsv,'largest','MaxIterations',50);
    
    % get new H
    H = U*S*V';
    
    % average antidiagonals of H to keep Hankel structure    
    H2 = fliplr(H);
    for ii=-(nrow-1):ncol-1
        vec = diag(H2,ii);
        avg = mean(vec);
        H2 = H2 + diag(avg - vec,ii);
    end
    H = fliplr(H2); % flip back
end
        
% get back fid
fid(1:ncol) = H(1,:);
fid(ncol:end) = H(:,end);
if doaug, fid = fid(1:end-1); end % remove extra point 
end
