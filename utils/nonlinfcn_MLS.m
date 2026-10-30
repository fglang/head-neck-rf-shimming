function [c, ceq, gc, gceq] = nonlinfcn_MLS(x, VOPs, SARLIM, stepdur, TR)
%NONLINFCN_MLS Evaluate VOP local-SAR constraints and their gradients.
%   [C, CEQ, GC, GCEQ] = NONLINFCN_MLS(X, VOPS, SARLIM, STEPDUR, TR)
%   supplies FMINCON inequality constraints for real/imaginary RF weights
%   X. VOPS contains local-SAR matrices, SARLIM is in W/kg, and STEPDUR and
%   TR are in seconds. There are no equality constraints.

arguments
    x
    VOPs
    SARLIM
    stepdur = 0.2e-3;
    TR = 10e-3;
end

% Explicit output branches let FMINCON request only the values it needs.
c = constraint_SAR(x, VOPs, SARLIM, stepdur, TR, 0);
if nargout > 1
    ceq = [];
end
if nargout > 2
    gc = constraint_SAR(x, VOPs, SARLIM, stepdur, TR, 1);
end
if nargout > 3
    gceq = [];
end
