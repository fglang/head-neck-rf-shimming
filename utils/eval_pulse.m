function [FAmap, pnrmse, ppower, SAR_VOP] = eval_pulse(p, A, mask, targFA, VOPs, swmask, stepdur, TR)
%EVAL_PULSE Evaluate an RF pulse in flip-angle, power, and SAR metrics.
%   [FAMAP, PNRMSE, PPOWER, SAR_VOP] = EVAL_PULSE(P, A, MASK, TARGFA,
%   VOPS, SWMASK, STEPDUR, TR) maps complex RF weights P through system
%   matrix A. MASK restores spatial dimensions; TARGFA and FAMAP are in
%   degrees; VOPS contains local-SAR matrices; SWMASK is retained for
%   compatibility with time-division code but is not used by this wrapper;
%   and STEPDUR and TR are in seconds. PNRMSE is fractional, PPOWER is
%   average power (W), and SAR_VOP is peak local SAR (W/kg).

arguments
    p % Complex pulse vector, N_t x N_c.
    A % Spatial-domain system matrix.
    mask % Spatial mask, N_x x N_y x N_z.
    targFA % Target flip angle [degrees].
    VOPs % VOP SAR matrices, N_c x N_c x N_VOP.
    swmask = [1] % TDM switch states: 1 for top row, 2 for bottom row.
    stepdur = 0.2e-3 % Time-step duration [s].
    TR = 10e-3; % Repetition time [s].
end

nc = size(VOPs,1);

FAmap = toSpatial(calcFA(p, A), mask);
pnrmse = nrmse2(targFA, abs(FAmap));
ppower = calcPowerForPulse(p, nc, stepdur, TR);

% Peak local SAR across the VOP set.
[~, SAR_VOP] = calcSarForPulse(p, VOPs, stepdur, TR);

end
