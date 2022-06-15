function [tau_0,t_a,MSE] = ...
    Dieterich_nucleation_sgolay_filter_test(order_range,npoints_range,...
        tau_0_range,t_a_range,chi2_delta,dt,rate_time,rate,geodetic_time,geodetic_V,n_pts_used)

tau_0 = {};
t_a = {};

for ii=1:length(order_range)
    for jj=1:length(npoints_range)
        disp([num2str(ii),' and ',num2str(jj)])
        
        geodetic_filtered = sgolayfilt(geodetic_V,order_range(ii),npoints_range(jj));
        geodetic_interpl = interp1(geodetic_time,geodetic_filtered,rate_time);
        
        [tau_0_out,t_a_out,mse_out] = ...
            Dieterich_nucleation_grid_search(tau_0_range,t_a_range,chi2_delta,dt,...
            rate,geodetic_interpl,0,n_pts_used);
        tau_0(ii,jj).val = tau_0_out.val;
%         tau_0(ii,jj).min = tau_0_out.min;
%         tau_0(ii,jj).max = tau_0_out.max;
        t_a(ii,jj).val = t_a_out.val;
%         t_a(ii,jj).min = t_a_out.min;
%         t_a(ii,jj).max = t_a_out.max;
        MSE(ii,jj) = mse_out;
    end
end
    
    
end