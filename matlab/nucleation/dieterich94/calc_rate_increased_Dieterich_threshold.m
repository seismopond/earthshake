function [time,RateIncrease] = calc_rate_increased_Dieterich_threshold(tau_0,tau_c,t_a,time,geodetic,dt)

%%% Initial creation: August 18, 2021

% Calculte Seismicity rate prediction for the best value 
for ii=1:length(time)
    Kernel(ii)=exp(tau_0*geodetic(ii)-tau_c);
    if ii==1
        Kernel_Integral(ii)=Kernel(ii)*(dt);
    else
        Kernel_Integral(ii)=Kernel(ii)*dt+Kernel_Integral(ii-1);
    end
    RateIncrease(ii)=Kernel(ii)/(1+(1/t_a)*Kernel_Integral(ii));
end

end