function res = mdl2_single_stochastic_sim(S,r, K_H,tmax)
% mdl1_single_stochastic_sim simulates stochastic trajectory for type I TA
% model 2
%   A single stochastic trajectory  is simulated from t = 0  to t = tmax, 
%   for given values of S, r and K_H and preset values of other model
%   parameters (see below)

    % Parameters and starting conditions:
    ht0 = 5;
    H0= 0;
    lambda = 1/3;
    lambda_d = 1/30; 
    kt = 10;
    n = 4;
    t0 = 0;

    % Initiate the time and molecule count vectors:
    t_vec = [];
    ht_vec = [];
    H_vec = [];
    t = t0;
    ht = ht0;
    H = H0;

    % Start the simulation:
                    
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
        p3 = kt*ht/(H/K_H +1) - kt*S*ht';
        p4 = H*lambda_d/((H/K_H)^n +1);
        propensities = [p1, p2, p3, p4];
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

    res = [t_vec; H_vec];

end