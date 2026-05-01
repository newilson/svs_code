function [fidecc, fid_refecc, phcorr] = eddyCurrentCorrection(fid_ref,fid,opt,nc)
%
% [fidecc, fid_refecc, phcorr] = eddyCurrentCorrection(fid_ref,fid,opt,nc)
%
% fid_ref is [npts x nc x other stuff]
% fid is [npts x nc x other stuff]
% opt: -1 BROWN, 0 ECC, 1 QUALITY, 2 QUECC
% nc: (optional) number of coils. If nc==1, inserts singleton coil
%     dimension so that dim 2 is not mistaken for averages.

if nargin<3 || isempty(opt)
    opt = 0;
end

% handle single-coil case where coil dimension was squeezed out
singleCoil = nargin >= 4 && nc == 1 && size(fid,2) ~= 1;
if singleCoil
    fid = reshape(fid, size(fid,1), 1, []);
    fid_ref = reshape(fid_ref, size(fid_ref,1), 1, []);
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
    case 0 % ECC: Klose MRM 1990 - with signal range restriction
        si = size(fid);
        npts = si(1);
        nc = si(2);
        fid = reshape(fid,npts,nc,[]);
        si_ref = size(fid_ref);
        fid_ref = reshape(fid_ref,npts,nc,[]);
        phcorr = exp(-1i*angle(fid_ref));

        % per-FID noise threshold and phase replacement (one column = one coil/channel FID)
        F2 = reshape(fid_ref, npts, []);
        P2 = reshape(phcorr,  npts, []);
        tailStart = max(npts-30, 1);
        noiseStd = std(F2(tailStart:end, :), 0, 1);
        for jj = 1:size(F2,2)
            noiseinds = find(abs(F2(:,jj)) < 5*noiseStd(jj));
            if ~isempty(noiseinds)
                P2(noiseinds, jj) = P2(max(noiseinds(1)-1,1), jj);
            end
        end
        phcorr = reshape(P2, npts, nc, []);

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

% remove singleton coil dimension if we inserted it
if singleCoil
    fidecc = squeeze(fidecc);
    fid_refecc = squeeze(fid_refecc);
    phcorr = squeeze(phcorr);
end

            