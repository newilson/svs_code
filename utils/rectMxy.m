function m = rectMxy(thDeg,Tus,dfHz,rel)
% RECTMXY  Transverse magnetization after an off-resonance rectangular pulse.
%
%   m = rectMxy(thDeg, Tus, dfHz, rel)
%
% M0 = z before the pulse; rotation by w_eff*T about the tilted effective
% field.  Exact Bloch solution, not the small-angle sinc approximation.
%
%   thDeg  prescribed (nominal, on-resonance) flip angle [deg]
%   Tus    pulse duration [us]
%   dfHz   offset of the observed peak from the excitation center [Hz]
%   rel    relB1 = achieved / prescribed B1 at that location
%
% Two distinct effects are folded in: rel is the SPATIAL transmit scaling at
% the location of the peak, dfHz the SPECTRAL excitation profile.  A hard
% pulse of duration T has its first profile null at |dfHz| = 1/T, so a long
% pulse (used to reach flips well past 90 deg) buys flip range at the cost of
% a narrow profile.
%
% Shared by isis_svs_refVolt (fitting relB1 from a flip series) and
% extRefQuant31P (applying the excitation factor at the acquisition's own
% pulse, center and flip).  Those two differ, so never reuse one's factor for
% the other -- recompute with the geometry that applies.

T  = Tus*1e-6;
w1 = rel*deg2rad(thDeg(:))/T;
dw = 2*pi*dfHz;
we = hypot(w1,dw);
we(we==0) = eps;
al = we*T;
sb = w1./we;
cb = dw./we;
m  = sb .* sqrt( cb.^2 .* (1-cos(al)).^2 + sin(al).^2 );
m  = reshape(m,size(thDeg));
end
