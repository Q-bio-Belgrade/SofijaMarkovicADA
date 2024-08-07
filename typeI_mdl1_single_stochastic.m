% Type I TA model 1: drawing single srochastic sim + low and high
% determnnistic state
% Author: Sofija Markovic
% Date: 08/12/2023

clear
%% parameters for stochastic simulations:
S_vec = 0.25;
r_vec = 0.5;
KH_vec = 100;
parfor i = 1:10
    for S=S_vec
        for r = r_vec
            for K_H = KH_vec
                % Starting conditions:
                ht0 = 5;
                H0_low = 0;
                H0_high = 200;
                lambda = 1/3;
                lambda_d = 1/30; 
                kt = 10;
                n = 4;
                tmax = 10000;
                t0 = 0;
                t_vec = [];
                ht_vec = [];
                H_vec = [];
                t = t0;
                ht = ht0;
                H = H0_low;
                propensities_vec = [];
                            
                while t < tmax 
                    
                    t_vec = [t_vec, t];
                    ht_vec = [ht_vec, ht];
                    if H <0
                        H = 0;
                    end
                    H_vec = [H_vec, H];
                
                    % Calculate the reaction propensities:
                    p1 = r*lambda*lambda_d*K_H/kt;
                    p2 = lambda*ht;
                    p3 = kt*ht - kt*S*ht;
                    p4 = kt*S*ht*H/K_H + H*lambda_d/((H/K_H)^n +1);
                    propensities = [p1, p2, p3, p4];
                    propensities_vec = [propensities_vec; propensities];
                    p_tot = sum(propensities);
                
                    % Get the probabilities for each of the reactions:
                    probs = propensities./p_tot;
                
                    % Probability distribution:
                    my_dist = makedist('Multinomial','Probabilities',probs);
                
                    % Waiting times:
                    dt = exprnd(1/p_tot);
                    t = t + dt;
                
                    switch random(my_dist)
                    case 1
                        ht = ht + 1 ;
                    case 2
                        ht = ht - 1;
                    case 3
                        H = H + 1;
                    case 4
                        H = H - 1;
                    end
                
                end
                
                plot(t_vec, H_vec,"LineWidth",1.5)
                
                hold on
                
                [T,Z] = ode45(@(t,z)rescaled_mdl1_odes_1(t,z,r,S,K_H, lambda_d),0:0.1:tmax,H0_low);
                
                plot(T,Z(:,1),"LineWidth",2)
                
                [T,Z] = ode45(@(t,z)rescaled_mdl1_odes_1(t,z,r,S,K_H, lambda_d),0:0.1:tmax,H0_high);
                
                plot(T,Z(:,1),"LineWidth",2)
                
                title(['r = ',num2str(r), ' S = ', num2str(S),' K_H =', num2str(K_H)])
                
                legend('stochastic H_0 = 0','deterministoc H_0 = 0', 'deterministic H_0 = 200')
                
                hold off
                fig_title = "TypeI_stochastic_r0"+ num2str(r*10) + "_S0" + num2str(S*100) + "_KH" + num2str(K_H) + "_" + num2str(i);
                saveas(gcf, fig_title, "fig")
                saveas(gcf, fig_title, "png")
            end
        end
    end
end
