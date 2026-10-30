function [p_fmincon_constr_kopt, nrmse_fmincon_kopt, VOPSAR_fmincon_kopt] = ...
    calcLcurve(A, targFA, p0, VOPs, SARLIMs, mask, doWarmStart, algo, stepdur, TR, maxiter)
%CALCLCURVE Solve magnitude-least-squares RF shims over SAR limits.
%   [P, NRMSE, PEAKSAR] = CALCLCURVE(A, TARGFA, P0, VOPS, SARLIMS, MASK,
%   DOWARMSTART, ALGO, STEPDUR, TR, MAXITER) optimizes one RF weight vector
%   per local-SAR limit. A is the small-tip-angle system matrix; TARGFA is
%   in degrees; P0 and P contain complex channel weights; VOPS contains the
%   local-SAR matrices; SARLIMS and PEAKSAR are in W/kg; MASK maps target
%   voxels back to space; STEPDUR and TR are in seconds. NRMSE is fractional.

arguments
    A
    targFA
    p0
    VOPs
    SARLIMs
    mask
    doWarmStart = true
    algo = 'interior-point'
    stepdur = 0.2e-3;
    TR = 10e-3;
    maxiter = 500;
end

if isempty(p0)
    p0 = randn(size(A,2),1) + 1j.*randn(size(A,2),1);
end

pulse_opt_constr = complex2ReIm(p0);

p_fmincon_constr_kopt = zeros(length(p0), length(SARLIMs));
nrmse_fmincon_kopt = zeros(length(SARLIMs),1);
VOPSAR_fmincon_kopt = zeros(length(SARLIMs),1);

for jj = 1:length(SARLIMs)
    objfun = @(x) loss_MLS(x, A, deg2rad(targFA));
    
    if doWarmStart
        x0 = pulse_opt_constr; % warm start
    else
        x0 = complex2ReIm(p0);
    end
    lb = [];
    ub = [];

    nonlinfcn = @(x) nonlinfcn_MLS(x, VOPs, SARLIMs(jj), stepdur, TR);

    options = optimoptions("fmincon", ...
        'Algorithm', algo, ...
        "Display","none", ...
        "MaxIterations", maxiter, ...
        "MaxFunctionEvaluations", 5000, ...
        "SpecifyObjectiveGradient", true, ...
        "SpecifyConstraintGradient", true);

    pulse_opt_constr = fmincon(objfun, x0, [], [], [], [], lb, ub, nonlinfcn, options);

    p_fmincon_constr_kopt(:,jj) = reIm2Complex(pulse_opt_constr);
    FAmap_fmincon_kopt = toSpatial(calcFA(p_fmincon_constr_kopt(:,jj), A), mask);
    nrmse_fmincon_kopt(jj) = nrmse2(targFA, abs(FAmap_fmincon_kopt));

    if nargout > 2
        [~, VOPSAR_fmincon_kopt(jj)] = calcSarForPulse( ...
            p_fmincon_constr_kopt(:,jj), VOPs, stepdur, TR, [], 1);
    end

    fprintf('L-curve: DONE %d/%d\n', jj, length(SARLIMs));
end

end
