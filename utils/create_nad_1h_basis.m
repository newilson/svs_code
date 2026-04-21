%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create_nad_1h_basis.m
%
% Generate 1H NAD+ downfield LCModel basis set (.RAW + .IN files) at 7T
%
% Based on de Graaf et al., Magn Reson Med 78:828-835 (2017)
% "Detection of cerebral NAD+ in humans at 7T"
%
% Basis set includes:
%   NAD_H2, NAD_H4, NAD_H6  - NAD+ singlet resonances
%   NAA                      - N-acetyl aspartate amide singlet
%   S1-S8                    - 8 unspecified singlet signals
%
% Chemical shifts from Table/text in de Graaf et al.:
%   NAD+ H2  = 9.33 ppm
%   NAD+ H6  = 8.83 ppm
%   NAD+ H4  = 9.14 ppm
%   NAA amide = 7.85 ppm
%   S1-S8    = 7.27, 7.51, 7.93, 8.04, 8.16, 8.27, 8.39, 8.50 ppm
%
% Usage:
%   1. Edit the CONFIGURATION section below
%   2. Run in MATLAB: create_nad_1h_basis
%   3. Build final .BASIS file on LCModel machine:
%        makebasis < nad_1h_downfield_7T.IN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% ===== CONFIGURATION =================================================

% Output directory for .RAW and .IN files
outputdir = fullfile(fileparts(mfilename('fullpath')), '..', 'basis', 'nad_1h_downfield_7T');

% Basis file name (no extension)
basisfilename = 'nad_1h_downfield_7T';

% Field strength in Tesla
Bo = 7;

% Exponential line broadening in Hz (default: 2)
lwHzDefault = 2;

%% ===== ACQUISITION PARAMETERS ========================================

% 1H gyromagnetic ratio: 42.576 MHz/T
hzpppm = 42.576 * Bo;  % = 297.2 MHz at 7T

acqparam = struct();
acqparam.hzpppm         = hzpppm;
acqparam.nbpoints       = 2048;
acqparam.sw_hz          = 12000;          % 12 kHz spectral width (matching de Graaf)
acqparam.receiveroffset_ppm = 0;          % will be set per-metabolite in createspectrum_1h

fprintf('Acquisition parameters (7T 1H):\n');
fprintf('  hzpppm     = %.4f\n', acqparam.hzpppm);
fprintf('  nbpoints   = %d\n', acqparam.nbpoints);
fprintf('  sw_hz      = %d\n', acqparam.sw_hz);
fprintf('  deltat     = %.6e s\n', 1/acqparam.sw_hz);

%% ===== SET UP FILE STRUCTURE =========================================

if ~exist(outputdir, 'dir')
    mkdir(outputdir);
    fprintf('\nCreated output directory: %s\n', outputdir);
end

if outputdir(end) ~= filesep
    outputdir = [outputdir filesep];
end

%% ===== DEFINE METABOLITES ============================================
% All resonances are singlets for the 1H downfield NAD+ basis set.
% Chemical shifts from de Graaf et al. MRM 2017.

clear metabolite
metabolite = struct('name', {}, 'cs', {}, 'ampl', {}, 'phase', {});

% NAD+ H2 singlet at 9.33 ppm
metabolite(end+1).name  = 'NADH2';
metabolite(end).cs      = 9.33;
metabolite(end).ampl    = 1;
metabolite(end).phase   = 0;

% NAD+ H4 singlet at 9.14 ppm
metabolite(end+1).name  = 'NADH4';
metabolite(end).cs      = 9.14;
metabolite(end).ampl    = 1;
metabolite(end).phase   = 0;

% NAD+ H6 singlet at 8.83 ppm
metabolite(end+1).name  = 'NADH6';
metabolite(end).cs      = 8.83;
metabolite(end).ampl    = 1;
metabolite(end).phase   = 0;

% NAA amide singlet at 7.85 ppm
metabolite(end+1).name  = 'NAA';
metabolite(end).cs      = 7.85;
metabolite(end).ampl    = 1;
metabolite(end).phase   = 0;

% 8 unspecified singlets (S1-S8) from de Graaf et al.
singlet_ppm = [7.27, 7.51, 7.93, 8.04, 8.16, 8.27, 8.39, 8.50];
for kk = 1:length(singlet_ppm)
    metabolite(end+1).name  = sprintf('S%d', kk);
    metabolite(end).cs      = singlet_ppm(kk);
    metabolite(end).ampl    = 1;
    metabolite(end).phase   = 0;
end

nbmetabolites = length(metabolite);

% Check metabolite name lengths (LCModel limit: 6 chars)
idx = find(arrayfun(@(m) length(m.name), metabolite) > 6);
if ~isempty(idx)
    arrayfun(@(m) fprintf('  %s (%d chars)\n', m.name, length(m.name)), metabolite(idx));
    error('Metabolite name(s) exceed 6 characters (LCModel limit)');
end

fprintf('\nBasis set metabolites (%d total):\n', nbmetabolites);
for ii = 1:nbmetabolites
    fprintf('  %6s : %.2f ppm\n', metabolite(ii).name, metabolite(ii).cs);
end

%% ===== GENERATE BASIS SPECTRA ========================================

% Reference position for createspectrum (far from peaks of interest)
tms_ref = 0;

lb = lwHzDefault;

for ii = 1:nbmetabolites
    createspectrum_1h( ...
        [metabolite(ii).cs, tms_ref], ...
        [metabolite(ii).ampl, 0.1], ...
        [metabolite(ii).phase, 0], ...
        lb, 0, outputdir, metabolite(ii).name, acqparam);
end

%% ===== CREATE .IN FILE FOR MAKEBASIS =================================

infile = [outputdir basisfilename '.IN'];
fprintf('\nWriting %s\n', infile);

fileid = fopen(infile, 'w');

fprintf(fileid, ' $NMALL\n');
fprintf(fileid, ' HZPPPM=%g\n', acqparam.hzpppm);
fprintf(fileid, ' NUNFIL=%d\n', acqparam.nbpoints);
fprintf(fileid, ' DELTAT=%g\n', 1/acqparam.sw_hz);
fprintf(fileid, ' FILBAS=''%s.BASIS''\n', basisfilename);
fprintf(fileid, ' FILPS=''%s.PS''\n', basisfilename);
fprintf(fileid, ' AUTOSC=.FALSE.\n');
fprintf(fileid, ' AUTOPH=.FALSE.\n');
fprintf(fileid, ' PPMST=12.\n');
fprintf(fileid, ' PPMEND=5.\n');
fprintf(fileid, ' $END\n\n');

for ii = 1:nbmetabolites
    fprintf(fileid, ' $NMEACH\n');
    fprintf(fileid, ' FILRAW=''%s.RAW''\n', metabolite(ii).name);
    fprintf(fileid, ' METABO=''%s''\n', metabolite(ii).name);
    fprintf(fileid, ' DEGZER=0.\n');
    fprintf(fileid, ' DEGPPM=0.\n');
    fprintf(fileid, ' CONC=1.\n');
    fprintf(fileid, ' PPMAPP=1.,-1.\n');
    fprintf(fileid, ' PPMPK=%.2f\n', metabolite(ii).cs);
    fprintf(fileid, ' $END\n\n');
end

fclose(fileid);

%% ===== SUMMARY =======================================================

fprintf('\n===== GENERATION COMPLETE =====\n\n');

rawfiles = dir(fullfile(outputdir, '*.RAW'));
fprintf('Generated %d .RAW files in %s:\n', length(rawfiles), outputdir);
for ii = 1:length(rawfiles)
    fprintf('  %s\n', rawfiles(ii).name);
end

infiles = dir(fullfile(outputdir, '*.IN'));
fprintf('\nGenerated %d .IN file(s):\n', length(infiles));
for ii = 1:length(infiles)
    fprintf('  %s\n', infiles(ii).name);
end

fprintf('\nNext step: copy basis directory to LCModel machine and run:\n');
fprintf('  cd %s\n', outputdir);
fprintf('  makebasis < %s.IN\n', basisfilename);
