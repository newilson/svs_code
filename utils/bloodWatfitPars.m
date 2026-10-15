function watfit = bloodWatfitPars()
% BLOODWATFITPARS  Blood 1H water-reference fit model, shared by the
% single-TE water reference (processU01_Blood_1H -> bru2mat_svsNAD watfit)
% and the water T2 series (processU01_Blood_waterT2): two co-located complex
% Voigts (narrow core + broad pedestal), no baseline, no lineshape kernel.

watfit = struct();
watfit.method        = 'fit';
watfit.ppm_range     = [4 5.5];
watfit.ph_range      = [-60 60];
watfit.mode          = 5;
watfit.nComp         = 2;
watfit.twoStagePhase = false;
watfit.centerInit    = [4.70 4.70];
watfit.centerBounds  = [4.64 4.76;
                        4.58 4.82];
watfit.widthInit     = [0.03 0.10];
watfit.baselineOpt   = struct('enable', false, 'knotSpacing', 0.5, 'lambda', 1, ...
                              'lambdaAmpl', 0, 'TolFun', 1e-6, ...
                              'widthBounds', [0.015 0.06;
                                              0.06  0.45]);
watfit.lineShapeOpt  = struct('enable', false, 'areaPadFactor', 100);
end
