function [A, optmask] = get_Amat_individ(B_plus, mask, ix, maskmode)
%GET_AMAT_INDIVID Build a single-subject small-tip-angle system matrix.
%   [A, OPTMASK] = GET_AMAT_INDIVID(B_PLUS, MASK, IX, MASKMODE) extracts
%   subject IX from B_PLUS (T/V) and MASK. MASKMODE = 1 uses the subject's
%   anatomical target; other values use the full field of view. A maps
%   complex RF channel weights to flip angle, and OPTMASK contains ones in
%   evaluated voxels and NaNs elsewhere.

if nargin < 4
    maskmode = 1;
end
if maskmode == 1
    optmask = mask(:,:,ix);
else
    optmask = ones(size(mask(:,:,ix)));
end
optmask = logtonan(optmask);
stepSize = [2,2,3].*1e-3;
stepdur = 1e-3;
dist = createVoxelGrid(optmask, stepSize) .* optmask;
dist = reshape(dist(~isnan(dist)), [], 3);
B1m_UP = double(B_plus(:,:,:,ix) .* optmask);
s = flattenize(B1m_UP, 3);
A = createAmat(s, dist, [0, 0, 0], stepdur);
end
