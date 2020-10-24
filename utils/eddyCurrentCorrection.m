function [fidecc, fid_refecc] = eddyCurrentCorrection(fid_ref,fid,opt)
%
% [fidecc, fid_refecc] = eddyCurrentCorrection(fid_ref,fid,opt)
%
% fid_ref is [npts x nc]
% fid is [npts x nc x other stuff]
% opt: 0 ECC, 1 QUALITY, 2 QUECC

if nargin<3 || isempty(opt)
    opt = 0; 
end

switch opt
    case 0 % ECC: Klose MRM 1990
        phcorr = exp(-1i*angle(fid_ref));
        fid_refecc = fid_ref .* phcorr;
        si = size(fid);
        npts = si(1);
        nc = si(2);
        fid = reshape(fid,npts,nc,[]);
        phcorr = repmat(phcorr,[1 1 size(fid,3)]);
        fidecc = fid .* phcorr;
        fidecc = reshape(fidecc,si);
    case 1 % QUALITY: deGraaf et al MRM 1990
        warning('not yet implemented')
    case 2 % QUECC: Bartha et al MRM 2000
        warning('not yet implemented')
end

            