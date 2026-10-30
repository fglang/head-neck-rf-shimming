function [f,g] = loss_MLS(xri, A, b, lbd_power)
%LOSS_MLS Evaluate the magnitude-least-squares objective and gradient.
%   [F, G] = LOSS_MLS(XRI, A, B, LBD_POWER) evaluates
%   ||abs(A*x)-abs(B)||_2^2 plus optional RF-power regularization. XRI stores
%   the real and imaginary parts of complex channel weights x; A is the
%   system matrix; B is the target flip angle in radians; LBD_POWER is the
%   regularization weight. F is scalar and G is ordered like XRI.

if nargin < 4
    lbd_power = 0;
end
x = reIm2Complex(xri);
Ax = A*x;
y = abs(Ax) - abs(b);
pulse_power = sum(abs(x).^2./50);
f = sum(y.^2, 'all') + lbd_power .* pulse_power;

if nargout > 1 % Analytical gradient.
    temp1 = exp(1j.*angle(Ax));
    temp = A'*(y.*temp1);
    gRe = 2*real(temp);
    gIm = 2*imag(temp);
    g = [gRe; gIm] + 2.*lbd_power.*xri ./ 50;
end
end
