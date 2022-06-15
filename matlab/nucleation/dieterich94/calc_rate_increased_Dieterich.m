function [time,RateIncrease] = calc_rate_increased_Dieterich(tau_0,t_a,time,geodetic,dt)

%%% Initial creation: August 18, 2021

% Calculte Seismicity rate prediction for the best value 
for ii=1:length(time)
    Kernel(ii)=exp(tau_0*geodetic(ii));
    if ii==1
        Kernel_Integral(ii)=Kernel(ii)*(dt);
    else
        Kernel_Integral(ii)=Kernel(ii)*dt+Kernel_Integral(ii-1);
    end
    RateIncrease(ii)=Kernel(ii)/(1+(1/t_a)*Kernel_Integral(ii));
end

end