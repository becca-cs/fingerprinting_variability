function [T] = output_numerical_hass3(N,c,lam,gam,alph,f)
    % Generates annual temperature timeseries for 3-box model
    % INPUTS 
    %   N   : number of years for simulation to run; longer simulations 
    %           (≥10,000 yrs) recommended for stability of power spectra
    %   c   : ocean heat capacity [W m^-2 K^-1 yr^-1]
    %   lam : climate feedback strength; should be positive [W m^-2 K^-1]
    %   gam : coupling parameter, determined by diffusivity [W m^-2 K^-1]
    %   alph: moisture constant. Setting > 1 allows for preferential
    %           transport of energy from tropics to extratropics. In
    %           current climatology, this value is ~1.7.
    %   f   : f(1) standard deviation of radiative anomalies in 
    %           extratropical Southern Hemisphere box 
    %         f(2) standard deviation of radiative anomalies in 
    %           tropical box 
    %         f(3) standard deviation of radiative anomalies in 
    %           extratropical Northern Hemisphere box 
    %         f(4) standard deviation of globally coherent radiative 
    %           anomalies
    %         Units are [W m^-2]
    % OUTPUTS
    %   T   : annual temperature [K] in each box. N x 3 array.

if lam<0
    disp('Lambda must be positive')
    return;
end

if numel(f)~=4
    disp('Number of elements in f must be 4.')
    return;
end

dt = 1/10; % timestep
F = f.*randn(N+1,4); % normally distributed forcing

% you can add an ENSO signal with period centered around 5 years
% comment out if don't want ENSO!
%     x = 0:N;
%     ENSO = 0*sin(2*pi*(1/5)*x);
%
%     F(:,2) = F(:,2) + ENSO';

% initialize
T = zeros(N+1,3);
T(1,:) = F(1,1:3);

% run forward N years
for i = 2:N+1

    temp = T(i-1,:);
    ct = 1;
    for j = 1:(1/dt)

        % heat exchanges
        H21 = gam*(alph*temp(2)-temp(1));
        H23 = gam*(alph*temp(2)-temp(3));

        % temperature change
        dT1 = 1/c*(-lam*temp(1)+H21+F(i,1)+F(i,4));
        dT2 = 1/c*(-lam*temp(2)-H21-H23+F(i,2)+F(i,4));
        dT3 = 1/c*(-lam*temp(3)+H23+F(i,3)+F(i,4));
        temp = temp+[dT1 dT2 dT3].*dt;
        ct = ct+1;
    end

    % update temperature
    T(i,:) = temp;
end

T = T(2:end,:);

end

