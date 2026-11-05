function Phi = rbf_kernel(R, kernel, s)
switch kernel
    case 'gaussian',        Phi = exp(-(R/s).^2);
    case 'multiquadric',    Phi = sqrt(1 + (R/s).^2);
    case 'invmultiquadric', Phi = 1 ./ sqrt(1 + (R/s).^2);
    case 'cubic',           Phi = R.^3;
    case 'thinplate',       Phi = R.^2 .* log(R + (R==0));
    otherwise, error('unknown kernel %s', kernel);
end
end
