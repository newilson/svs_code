function appendPngsToPptx(pptxFile, pngFiles, titles)
% APPENDPNGSTOPPTX  Append PNG figures to a PowerPoint deck.
%
%   appendPngsToPptx(pptxFile, pngFiles)
%   appendPngsToPptx(pptxFile, pngFiles, titles)
%
% pngFiles is a cell array of paths or one wildcard pattern, e.g.
% 'C:\Users\CAMIPM-NW\Documents\859875_U01\SummaryResults\S01\S01_*.png'.
%
% Files named like the compareSubjectResultsU01 outputs
%   <sid>_<anatomy>_table.png
%   <sid>_<anatomy>_v<k>_<yyyy-mm-dd>_1H_fine_fit_region<n>_<name>.png
%   <sid>_<anatomy>_v<k>_<yyyy-mm-dd>_31P_fit.png
%   <sid>_<anatomy>.png                     (overlay)
% are ordered by subject, then anatomy (brain, calf, blood), then table, fit
% figures (1H regions in order, then 31P) and overlay.  The same fit figure
% from different visits shares one slide, side by side in visit order, each
% captioned with its visit and date (a fit figure from a single visit still
% gets its caption).  Every other file gets a slide of its
% own after these, in the order given (alphabetical for a pattern).
%
% Each slide uses the "Title Only" layout (the deck must have it), with the
% title box shrunk to its text and the pictures scaled to fit below it.
% titles (cell array, one per slide) defaults to the file name with
% underscores as spaces, without the visit and date on a shared slide.
%
% The deck is edited by PowerPoint itself (COM, Windows only), so its theme,
% layouts and existing slides are untouched.  If the deck is already open in
% PowerPoint the slides are added to the open copy and it is saved.  If it
% does not exist it is created (16:9, PowerPoint's default theme), along
% with its folder.

if nargin < 3, titles = {}; end
pptxFile = absPath(pptxFile);
if ~isfolder(fileparts(pptxFile)), mkdir(fileparts(pptxFile)); end
if ischar(pngFiles) || (isstring(pngFiles) && isscalar(pngFiles))
    d = dir(char(pngFiles));
    pngFiles = fullfile({d.folder}, {d.name});
end
pngFiles = cellfun(@absPath, cellstr(pngFiles), 'UniformOutput', false);
assert(~isempty(pngFiles), 'no PNG files given');
assert(all(cellfun(@isfile, pngFiles)), 'PNG not found: %s', strjoin(pngFiles(~cellfun(@isfile, pngFiles)), ', '));
[groups, autoTitles, captions] = slidePlan(pngFiles);
if isempty(titles), titles = autoTitles; end
titles = cellstr(titles);
assert(numel(titles) == numel(groups), 'need one title per slide (%d)', numel(groups));

app = actxserver('PowerPoint.Application');
pres = [];
for i = 1:app.Presentations.Count
    p = app.Presentations.Item(i);
    if strcmpi(p.FullName, pptxFile), pres = p; end
end
if isempty(pres) && ~isfile(pptxFile)
    pres = app.Presentations.Add(0);
    cleanup = onCleanup(@() closeDeck(app, pres)); %#ok<NASGU>
    pres.PageSetup.SlideWidth  = 960;   % 13.33 x 7.5 in (16:9)
    pres.PageSetup.SlideHeight = 540;
    pres.SaveAs(pptxFile, 24);          % ppSaveAsOpenXMLPresentation
    fprintf('Created %s\n', pptxFile);
elseif isempty(pres)
    pres = app.Presentations.Open(pptxFile, 0, 0, 0);
    cleanup = onCleanup(@() closeDeck(app, pres)); %#ok<NASGU>
end

layout = [];
L = pres.SlideMaster.CustomLayouts;
names = cell(1, L.Count);
for i = 1:L.Count
    names{i} = L.Item(i).Name;
    if strcmp(names{i}, 'Title Only'), layout = L.Item(i); end
end
assert(~isempty(layout), 'no "Title Only" layout in %s (layouts: %s)', pptxFile, strjoin(names, ', '));

W = pres.PageSetup.SlideWidth;
H = pres.PageSetup.SlideHeight;
m = 0.3 * 72;
for k = 1:numel(groups)
    s = pres.Slides.AddSlide(pres.Slides.Count + 1, layout);
    t = s.Shapes.Title;
    t.TextFrame.TextRange.Text = titles{k};
    t.TextFrame.AutoSize = 1;   % ppAutoSizeShapeToFitText
    top = t.Top + t.Height + m/2;
    n = numel(groups{k});
    capH = 0;
    if any(~cellfun(@isempty, captions{k})), capH = 24; end
    cellW = (W - (n + 1)*m) / n;
    for j = 1:n
        f = groups{k}{j};
        x0 = m + (j - 1)*(cellW + m);
        if capH > 0
            tb = s.Shapes.AddTextbox(1, x0, top, cellW, capH);   % msoTextOrientationHorizontal
            tb.TextFrame.TextRange.Text = captions{k}{j};
            tb.TextFrame.TextRange.Font.Size = 16;
            tb.TextFrame.TextRange.ParagraphFormat.Alignment = 2;   % ppAlignCenter
        end
        info = imfinfo(f);
        sc = min(cellW / info.Width, (H - top - capH - m) / info.Height);
        w = info.Width * sc;
        h = info.Height * sc;
        pic = s.Shapes.AddPicture(f, 0, -1, x0 + (cellW - w)/2, top + capH, w, h);
        pic.AlternativeText = strtrim([titles{k} ' ' captions{k}{j}]);
    end
end
pres.Save;
fprintf('Added %d slide(s) to %s\n', numel(groups), pptxFile);
end

function [groups, titles, captions] = slidePlan(files)
[~, names] = cellfun(@fileparts, files, 'UniformOutput', false);
anatomies = {'brain', 'calf', 'blood'};
n = numel(files);
order = repmat([Inf 0 0 0 0], n, 1);
order(:,6) = (1:n)';
key = names(:)';
title1 = strrep(names(:)', '_', ' ');
cap = repmat({''}, 1, n);
for i = 1:n
    tk = regexp(names{i}, '^S(\d+)_(brain|calf|blood)((?:_.+)?)$', 'tokens', 'once');
    if isempty(tk), continue; end
    head = [str2double(tk{1}) find(strcmp(anatomies, tk{2}))];
    if isempty(tk{3})
        order(i,1:5) = [head 3 0 0];
    elseif strcmp(tk{3}, '_table')
        order(i,1:5) = [head 1 0 0];
    else
        fit = regexp(tk{3}, '^_v(\d+)_(\d{4}-\d{2}-\d{2})_(1H_fine_fit_region\d+_\w+|31P_fit)$', 'tokens', 'once');
        if isempty(fit), continue; end
        r = regexp(fit{3}, '^1H_fine_fit_region(\d+)_', 'tokens', 'once');
        sub = 1e6;
        if ~isempty(r), sub = str2double(r{1}); end
        order(i,1:5) = [head 2 sub str2double(fit{1})];
        key{i} = sprintf('S%s_%s_%s', tk{1}, tk{2}, fit{3});
        title1{i} = strrep(key{i}, '_', ' ');
        cap{i} = sprintf('Visit %s (%s)', fit{1}, fit{2});
    end
end
[~, o] = sortrows(order);
groups = {}; titles = {}; captions = {};
for i = o(:)'
    if ~isempty(groups) && order(i,3) == 2 && strcmp(key{i}, lastKey)
        groups{end}{end+1} = files{i};
        captions{end}{end+1} = cap{i};
    else
        groups{end+1} = files(i); %#ok<AGROW>
        titles{end+1} = title1{i}; %#ok<AGROW>
        captions{end+1} = cap(i); %#ok<AGROW>
    end
    lastKey = key{i};
end
end

function closeDeck(app, pres)
pres.Close;
if app.Presentations.Count == 0, app.Quit; end
end

function p = absPath(p)
p = char(p);
if ~(numel(p) > 1 && (p(2) == ':' || startsWith(p, '\\')))
    p = fullfile(pwd, p);
end
end
