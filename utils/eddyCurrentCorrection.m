function [fidecc, fid_refecc] = eddyCurrentCorrection(fid_ref,fid,opt)
%
% [fidecc, fid_refecc] = eddyCurrentCorrection(fid_ref,fid,opt)
%
% fid_ref is [npts x nc x other stuff]
% fid is [npts x nc x other stuff]
% opt: 0 ECC, 1 QUALITY, 2 QUECC

if nargin<3 || isempty(opt)
    opt = 0; 
end

switch opt
    case -1 % no ECC: Brown MRM
        si = size(fid);
        npts = si(1);
        nc = si(2);
        fid = reshape(fid,npts,nc,[]);
        si_ref = size(fid_ref);
        fid_ref = reshape(fid_ref,npts,nc,[]);
        phcorr = exp(-1i*squeeze(angle(fid_ref(1,:,:))));
        if length(si_ref)>2
            phcorr = permute(phcorr,[3 1 2]);
            fid_refecc = fid_ref .* repmat(phcorr,[npts 1 1]);
            fid_refecc = reshape(fid_refecc,si_ref);
            phcorr = repmat(phcorr,[npts 1 1]);
        else
            fid_refecc = fid_ref .* repmat(phcorr,[npts,1]);
            phcorr = repmat(phcorr,[npts 1 size(fid,3)]);
        end
        fidecc = fid .* phcorr;
        fidecc = reshape(fidecc,si);      
    case 0 % ECC: Klose MRM 1990
        si = size(fid);
        npts = si(1);
        nc = si(2);
        fid = reshape(fid,npts,nc,[]);
        si_ref = size(fid_ref);
        fid_ref = reshape(fid_ref,npts,nc,[]);
        phcorr = exp(-1i*angle(fid_ref));
        fid_refecc = fid_ref .* phcorr;
        fid_refecc = reshape(fid_refecc,si_ref);
        if length(si_ref)<=2
            phcorr = repmat(phcorr,[1 1 size(fid,3)]);
        end
        fidecc = fid .* phcorr;
        fidecc = reshape(fidecc,si);
    case 1 % QUALITY: deGraaf et al MRM 1990
        warning('not yet implemented')
    case 2 % QUECC: Bartha et al MRM 2000
        warning('not yet implemented')
end

            