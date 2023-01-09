function[pout, Vout] = einschrk(pin, bound, Vin)

lo = bound.lo; hi = bound.hi; welche = bound.which;

if nargin < 3
    trans = sqrt((hi-pin) ./ (pin-lo)); pout = (1-welche) .* pin + welche .* trans;
    Vout = [ ];
else
    trans = (hi+lo .* pin .^ 2) ./ (1 + pin.^2); pout = (1 - welche) .* pin + welche .* trans;
    % now adjust the standard errors
    trans = 2 * pin .* (lo - hi) ./ (1 + pin.^2).^2;
    d = (1 - welche) + welche .* trans; % eitherunity or delta method .
    J = diag(d); Vout = J * Vin * J;
end