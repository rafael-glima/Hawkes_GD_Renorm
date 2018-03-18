function g = KernelFunc( t, weight, freq, s, pattern, decayr, p, q)

% compute sine function-based kernel
% t: time
% weight: amplitude
% freq: frequency
% s: time shift
% pattern: 'sine' and 'square'

if nargin < 8
	q = 0.5
end

switch pattern
    case 'sine'
        % g = weight.*exp(-decayr*t).*0.5.*(1 + cos(2*2*pi*freq*t+pi*s ));
        g = weight*sin(pi*freq*t+pi*s );
    case 'square'
        g = weight*ceil(0.5*(1 - cos(2*pi*freq*t+pi*s )));
    case 'exponential'
        g = weight*exp(-decayr*t);
    case 'powerlaw'
        %cte = ((1+t)^-p)
        %g = weight*((1+t)^-p);
        g = weight*power(1.+t,-p);
    case 'rayleigh'
        g = weight*exp(-decayr*power(t,2))
    case 'q-exponential'
        if (q ~= 1.)
            g = weight*(1+(1-q).*t).^(1/(1-q))
	    g(g<0.) == 0.
        elseif (q==1.)
            g = weight*exp(-t)
        else
            g = 0.	
        end
    case 'other'
        g = weight*(1 - cos(2*pi*freq*t+pi*s )).*exp(-decayr*t);
    otherwise
        disp('Please assign a kernel function!');
end
