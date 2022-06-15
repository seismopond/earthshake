function [tau_0,t_a,MSE] = ...
    Dieterich_nucleation_grid_search(tau_0_range,t_a_range,chi2_delta,dt,...
        rate,geodetic,plot_output,n_pts_used)
    
if ~exist('n_pts_used','var')
    n_pts_used = length(rate);
end

% Generate a grid
[X,Y] = meshgrid(tau_0_range,t_a_range);
Z = zeros(length(t_a_range),length(tau_0_range));

% Do a grid search
for xx=1:length(t_a_range)
    for yy=1:length(tau_0_range)
        for ii=1:length(rate)
            Kernel(ii)=exp(X(xx,yy)*geodetic(ii));
            if ii==1
                Kernel_Integral(ii)=Kernel(ii)*(dt);
            else
                Kernel_Integral(ii)=Kernel(ii)*dt+Kernel_Integral(ii-1);
            end
            RateIncrease(ii)=Kernel(ii)/(1+(1/Y(xx,yy))*Kernel_Integral(ii));          
        end
        RateIncrease = RateIncrease(:);
        rate = rate(:);
        Z(xx,yy) = mean((RateIncrease(1:n_pts_used) - rate(1:n_pts_used)).^2);
    end
end

% Find best fit value (minimum mean squared)
MSE = min(Z(:));
[x,y]=find(Z==MSE);
tau_0.val = X(x,y);
t_a.val = Y(x,y);

% Calculate chi-squared (normalized by the best fit)
chi2 = Z/MSE;

% Find error ellipse
% tau_0.min = interp1(chi2(x,1:y),tau_0_range(1:y),1+chi2_delta);
% tau_0.max = interp1(chi2(x,y:end),tau_0_range(y:end),1+chi2_delta);
% t_a.min = interp1(chi2(1:x,y),t_a_range(1:x),1+chi2_delta);
% t_a.max = interp1(chi2(x:end,y),t_a_range(x:end),1+chi2_delta);


% Plot grid search output, if the user wants
if plot_output == 1
    figure;
    contourf(X,Y,chi2,50); hold on;
    h = colorbar;
    colormap(flipud(parula))
    ylabel(h,'Reduced Chi-squared','FontSize',12)
    xlabel('$\tau_0/(a\sigma)$','Interpreter','latex')
    ylabel('$t_a$ (hours)','Interpreter','latex')
    scatter(tau_0.val,t_a.val,100,'xk')
    set(gca,'Fontsize',12)
end


end