function [refVolt,output] = runIsisRefVolt(acq,pars)
%
% [refVolt,output] = runIsisRefVolt(acq,pars)
%
% Wrapper around isis_svs_refVolt for the 859875_S01_09012026 brain session.
% acq picks one of the two ISIS_SVS flip angle series in E100.  They differ
% only in the excitation center frequency, so PCr must come out the same from
% both; that agreement is the check on the fit.
%
%   1 or 'pcr'   series 0006-0010, excitation centerd at -4.7 ppm.  PCr sits
%                565 Hz off center, so the rectangular pulse profile matters
%                and the flip angle curve rises to 130 deg instead of peaking.
%                The external reference at 8 ppm is 1528 Hz off center, near
%                the first null of the 600 us pulse, and is not usable here.
%
%   2 or 'ref'   series 0011-0015, excitation centerd at 0 ppm.  PCr is on
%                resonance and the external sodium diphosphate reference at
%                8.18 ppm is 984 Hz off center.  refVolt is taken from PCr.
%
% Fields set in pars override the per acquisition defaults below.
%
% 09/2026 - NW

if nargin<1 || isempty(acq), acq = 1; end
if nargin<2 || isempty(pars), pars = struct(); end

dataDir = ['Z:\Data\Spectroscopy\U01_temp\859875_S01_09012026_brain' filesep ...
           '859875_S01_09012026_brain' filesep 'E100'];

switch lower(string(acq))
    case {"1","pcr"}
        name       = 'PCr only, excitation center -4.7 ppm';
        series     = 6:10;
        def.peaks  = { 'PCr', 0.00, 0.6 };
        def.combine= 'PCr';
    case {"2","ref","extref"}
        name       = 'PCr + external reference, excitation center 0 ppm';
        series     = 11:15;
        def.peaks  = { 'PCr', 0.00, 0.6
                       'Ref', 8.18, 0.5 };
        def.combine= 'PCr';
    otherwise
        error('acq must be 1 (''pcr'') or 2 (''ref''), not %s',string(acq))
end

% shared across both acquisitions; extraction mode and dropPts come from the
% isis_svs_refVolt defaults
def.plt = true;

pars = mergePars(def,pars);

files = findSeries(dataDir,series);

fprintf('\n================================================================\n');
fprintf('acquisition %s: %s\n',string(acq),name);
fprintf('series %s\n',mat2str(series));
fprintf('================================================================\n');

[refVolt,output] = isis_svs_refVolt(files,pars);
output.acq     = acq;
output.acqName = name;
output.series  = series;

end

% =====================================================================
function files = findSeries(dataDir,series)
% One spectroscopy file per series number, matched on the trailing 4 digits of
% the series folder name.
assert(isfolder(dataDir),'data folder not found: %s',dataDir);
files = cell(1,numel(series));
for ii=1:numel(series)
    d = dir(fullfile(dataDir,sprintf('*_%04d',series(ii))));
    d = d([d.isdir]);
    assert(~isempty(d),'no folder ending in _%04d under %s',series(ii),dataDir);
    assert(numel(d)==1,'%d folders end in _%04d under %s',numel(d),series(ii),dataDir);
    f = dir(fullfile(d.folder,d.name,'*.*'));
    f = f(~[f.isdir]);
    keep = ~cellfun('isempty',regexpi({f.name},'\.(ima|dcm|dat)$','once'));
    f = f(keep);
    assert(~isempty(f),'no .IMA/.dcm/.dat in %s',fullfile(d.folder,d.name));
    assert(numel(f)==1,'%d spectroscopy files in %s',numel(f),fullfile(d.folder,d.name));
    files{ii} = fullfile(f.folder,f.name);
end
end

% =====================================================================
function p = mergePars(def,usr)
% usr wins, one level of nesting deep.
p = def;
for f = fieldnames(usr)'
    if isstruct(usr.(f{1})) && isfield(def,f{1}) && isstruct(def.(f{1}))
        p.(f{1}) = mergePars(def.(f{1}),usr.(f{1}));
    else
        p.(f{1}) = usr.(f{1});
    end
end
end
