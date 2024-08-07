function hdot = rescaled_mdl1_odes_1(t,h,r,S,KH,lambda_d)
hdot = r.*KH.*lambda_d.*(1-S.*(1+h(1)./KH)) - lambda_d.*h(1)./(1 + (h(1)./KH)^4);
end