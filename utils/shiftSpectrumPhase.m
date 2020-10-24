function spec_ph = shiftSpectrumPhase(spec,ph,varargin)
%
% spec_ph = shiftSpectrumPhase(spec,ph);
%
% spec is [npts x whatever]
% ph is either a constant phase, 2 elements ([pc0 pc1]), [npts x 1], or the same size as spec

si = size(spec);
lph = length(ph);
if lph==1
    spec_ph = spec * exp(-1i*ph);
elseif lph==2 && length(varargin)==2
    pc0 = ph(1);
    pc1 = ph(2);
    % additional inputs
    ppm = varargin{1};
    pivot = varargin{2};
    
    linph = 2*pc1/abs(ppm(end)-ppm(1)) * (ppm-pivot);
    phasevec = exp(-1i*pi/180*(pc0 + linph));
    phasearr = repmat(phasevec,[1,si(2)]);
    spec_ph = spec .* phasearr;
elseif isvector(ph)
    if ~isequal(si(1),lph)
        error('inconsistent dimensions')
    end
    spec = reshape(spec,si(1),[]);
    ph = repmat(ph,[1 si(2)]);
    spec_ph = spec .* exp(-1i*ph);
    spec_ph = reshape(spec_ph,si);
elseif size(ph)==si
    spec_ph = spec .* exp(-1i*ph);
else
    error('inconsistent dimensions')
end
