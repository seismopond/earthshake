function seis_proj = seis_projection(seis,fault_model)

%%%% Aug 10, 2021: initial creation for strike-slip fault only

theta = 90-fault_model.strike(1);
m1 = tand(theta);
m2 = -1/m1;
c1 = fault_model.y(1) - m1*fault_model.x(1);
for ii=1:length(seis.x)
    c2 = seis.y(ii) - m2*seis.x(ii);
    seis_proj.x(ii) = (c2-c1)/(m1-m2); 
    seis_proj.y(ii) = m1*seis_proj.x(ii) + c1;
    seis_proj.z(ii) = seis.dep(ii);
end

seis_proj.x = seis_proj.x(:);
seis_proj.y = seis_proj.y(:);
seis_proj.z = seis_proj.z(:);

end