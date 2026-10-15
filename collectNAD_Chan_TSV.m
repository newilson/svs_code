function T = collectNAD_Chan_TSV(dataRoot, outPath)
% COLLECTNAD_CHAN_TSV  Gather NAD_Chan 1H results from every session into one TSV.
%
%   collectNAD_Chan_TSV(dataRoot)
%   collectNAD_Chan_TSV(dataRoot, outPath)
%
% Reads <session>\processed\<session>_<scan>_1H_abs.mat (scan = slab, vox)
% for every session folder under dataRoot, plus the matching brainseg.sh TSV
% (<seg>_svsmask_<scan>.tsv) and imscribe mask (svsmask_<scan>.nii).
% Default outPath: <dataRoot>\NAD_Chan_1H.tsv.  One row per session x scan.
%
% COLUMNS
%   session, subject, scan, TE_ms, TR_ms, nAvg
%   uM_NADH2, uM_NADH6, uM_NADH4, uM_Trp
%   f_CSF, f_GM, f_WM            fractions of CSF+GM+WM inside the mask
%   vol_CSF_mL, vol_GM_mL, vol_WM_mL, vol_GMWM_mL
%   vol_mask_mL                  whole mask, including non-brain voxels
%   SNR_NADH2                    fitted NADH2 peak height / std of the real
%                                spectrum over 11.5-13.5 ppm (linear trend
%                                removed), after the pipeline's lb_met
%   water_FWHM_Hz                fitted water-reference linewidth

arguments
    dataRoot (1,:) char
    outPath  (1,:) char = fullfile(dataRoot, 'NAD_Chan_1H.tsv')
end

scans = {'slab','vox'};
d = dir(dataRoot); d = d([d.isdir] & ~ismember({d.name},{'.','..'}));

rows = {};
for i = 1:numel(d)
    sess = d(i).name;
    pd   = fullfile(dataRoot, sess, 'processed');
    for j = 1:numel(scans)
        f = fullfile(pd, sprintf('%s_%s_1H_abs.mat', sess, scans{j}));
        if ~isfile(f), continue; end
        rows{end+1} = buildRow(sess, scans{j}, f, pd); %#ok<AGROW>
        fprintf('  %-28s %-4s\n', sess, scans{j});
    end
end
assert(~isempty(rows), 'No *_1H_abs.mat found under %s', dataRoot);

T = struct2table([rows{:}], 'AsArray', true);
T = sortrows(T, {'subject','scan'});
writetable(T, outPath, 'FileType','text', 'Delimiter','\t');
fprintf('Wrote %s  (%d rows)\n', outPath, height(T));
end

% =====================================================================
function r = buildRow(sess, scan, matFile, pd)
m = load(matFile, 'output', 'absQ', 'segFracs');

r.session = sess;
r.subject = subjectOf(sess);
r.scan    = scan;
r.TE_ms   = m.output.met.TE_ms;
r.TR_ms   = m.output.met.TR_ms;
r.nAvg    = m.output.met.nAvg;

for k = {'NADH2','NADH6','NADH4','Trp'}
    r.(['uM_' k{1}]) = NaN;
end
for k = 1:numel(m.absQ.metabolites)
    r.(['uM_' m.absQ.metabolites(k).name]) = m.absQ.metabolites(k).conc_mM * 1e3;
end

r.f_CSF = m.segFracs(1); r.f_GM = m.segFracs(2); r.f_WM = m.segFracs(3);

r.vol_CSF_mL = NaN; r.vol_GM_mL = NaN; r.vol_WM_mL = NaN; r.vol_GMWM_mL = NaN;
tsv = dir(fullfile(pd, sprintf('*_seg_svsmask_%s.tsv', scan)));
if ~isempty(tsv)
    S  = readtable(fullfile(tsv(1).folder, tsv(1).name), 'FileType','text', ...
                   'Delimiter','\t', 'VariableNamingRule','preserve');
    nm = string(S.labelname);
    v  = S.('volume(mm^3)') / 1e3;
    r.vol_CSF_mL  = sum(v(nm == "CSF"));
    r.vol_GM_mL   = sum(v(nm == "GM"));
    r.vol_WM_mL   = sum(v(nm == "WM"));
    r.vol_GMWM_mL = r.vol_GM_mL + r.vol_WM_mL;
end

r.vol_mask_mL = NaN;
mk = fullfile(pd, sprintf('svsmask_%s.nii', scan));
if isfile(mk)
    info = niftiinfo(mk);
    r.vol_mask_mL = nnz(niftiread(info) > 0) * prod(info.PixelDimensions(1:3)) / 1e3;
end

r.SNR_NADH2     = snrNADH2(m.output);
r.water_FWHM_Hz = waterFWHM(m.output);
end

% =====================================================================
function s = snrNADH2(output)
f   = output.met.fit;
ppm = output.met.ppm(:);
sp  = real(f.spec_full(:));
sel = ppm > 11.5 & ppm < 13.5;
x   = ppm(sel);
y   = sp(sel);
res = y - polyval(polyfit(x, y, 1), x);
iH2 = strcmp(f.names, 'NADH2');
s   = max(real(f.comp(:, iH2))) / std(res);
end

function w = waterFWHM(output)
wf = output.wat.fit;
fHz = (output.wat.hz(2) - output.wat.hz(1)) / (output.wat.ppm(2) - output.wat.ppm(1));
switch numel(wf.pars)
    case 6, w = wf.pars(6) * fHz;     % Voigt: amp pos widthL ph widthG width
    case 4, w = wf.pars(3) * fHz;     % Lorentzian: amp pos width ph
    otherwise, w = NaN;
end
end

function s = subjectOf(sess)
p = strsplit(sess, '_');
k = find(~cellfun('isempty', regexp(p, '^S\w+$', 'once')), 1);
if isempty(k), s = sess; else, s = p{k}; end
end
