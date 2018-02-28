function v = gkernel( t, time, landmark, sigma )
    dt = t-time;
    dis = repmat(dt, [length(landmark), 1]) - repmat(landmark(:), [1,length(time)]);
    Sigma = repmat(sigma(:), [1,length(time)]);
    %v = 1./(sqrt(2*pi).*Sigma) .* exp(-(dis.^2)./(2*Sigma.^2));
    v = exp(-(dis.^2)./(2*Sigma.^2));
end