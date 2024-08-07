function [H_low_ss,H_high_ss] = get_model2_steady_states(r,S,KH,lambda_d)
%get_model1_steady_states returns low and high steady state concentration of 
% toxin from model 1 type I TA system. 
%
% function solves model1 ODE for given parameter values of r, S and KH. It
% starts from H0 = 0, to het the low toxin steady state and than gradually
% increases H until the other (high toxin) steady state is obtained, or the
% loop breaking condition is achieved (to prevent infinite loop in case the
% the system is not bistable). 

H0 = 100;

[~, H_low] = ode45(@(t,h)rescaled_mdl2_odes_1(t,h,r,S,KH,lambda_d),0:10000, 0);
H_low_ss = H_low(end,1);

[~, H_high] = ode45(@(t,h)rescaled_mdl2_odes_1(t,h,r,S,KH,lambda_d),0:10000, H0);
H_high_ss = H_high(end,1);

round = 1;
while H_high_ss - H_low_ss < H_high_ss*0.1
    H0 = H0*2;
    [~, H_high] = ode45(@(t,h)rescaled_mdl2_odes_1(t,h,r,S,KH,lambda_d),0:10000, H0);
    H_high_ss = H_high(end,1);
    round = round + 1; 
    if round >10
        break
    end
end

display(['H0 = ' num2str(H0)])

end

