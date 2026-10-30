function out = logtonan(inp)
%LOGTONAN Convert a logical mask to double with zeros represented by NaN.
%   OUT = LOGTONAN(INP) preserves nonzero values and replaces zeros by NaN.
%   OUT is dimensionless and has the same dimensions as INP.

out = double(inp);
out(out == 0) = NaN;
end
