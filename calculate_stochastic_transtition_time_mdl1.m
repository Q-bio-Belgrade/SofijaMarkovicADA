%% Parameter tuning for stochastic simulations

% Sofija Markovic
% 11.04.2024. 

% Aim of the code:
% For a set of (r,S) combinations which result in bistable system,
% determine average transition times from low to high toxin state and back.

clear

S = 0.25;
r = 0.55;
K_H = 80;
lambda = 1/3;
lambda_d = 1/30; 
kt = 10;
n = 4;
tmax = 100000;
colormap = [77, 190, 238; 10, 140, 105; 162, 20, 47]/255;

res = mdl1_single_stochastic_sim(S,r,K_H,tmax);

figure
plot(res(1,:), res(2,:), 'Color',colormap(2,:))
hold on 

[T,Z_low] = ode45(@(t,z)rescaled_mdl1_odes_1(t,z,r,S,K_H, lambda_d),0:0.1:tmax,0);

plot(T,Z_low(:,1), 'Color',colormap(1,:))

[T,Z_high] = ode45(@(t,z)rescaled_mdl1_odes_1(t,z,r,S,K_H, lambda_d),0:0.1:tmax,300);

plot(T,Z_high(:,1), 'Color',colormap(3,:))

sim_num = num2str(randi(1000));

xlabel('t[min]')
ylabel('H')
filename = 'mdl1_long_stochastic_r06S020KH_80_'+ sim_num + '.mat';
save(filename, "res")

%% Smoothing filter

t = res(1,:);
x = res(2,:);
figure
yMedFilt = medfilt1(x,1500,'truncate');
plot(t,x, ...
     t,yMedFilt)
hold on
y_gaussian_smooth = smoothdata(x, "gaussian");
plot(t,y_gaussian_smooth)
legend('original signal','median filter', 'gaussian')
%% 
res(3,:) = y_gaussian_smooth; 

%% Distinguish high and low states in stochastic regime:

% Get low- and high-toxin steady states:

[H_low_ss, H_high_ss] = get_model1_steady_states(r,S,K_H,lambda_d);

% Choose whetrer to use the filter or not
row = 3;

for i = 1:size(res,2)    
    H = res(row,i);
    H_diff_low = abs(H-H_low_ss);
    H_diff_high = abs(H-H_high_ss);
    if H_diff_low < H_diff_high
        res(4,i) = 1;
    else
        res(4,i) = 2;
    end
end

%% Plot to check:
figure
gscatter(res(1,:), res(2,:),res(4,:))
legend('Low','High')

%% Calculate transition times:

% create empty cell for low and high state results:
low_cell = {};
high_cell = {};

i  = 1;
l = 1;
h = 1;

while i < size(res,2)
    low = [];
    high = [];

    while res(4,i) == 1
       low = [low, res(:,i)];
       i = i + 1;
       if i > size(res,2)
           break
       end
    end

    if isempty(low) == 0
        low_cell{1,l} = low;
        l = l + 1; 
    end

    if i > size(res,2)
        break
    end

    while res(4,i) == 2
        high = [high, res(:,i)];
        i = i+1;
        if i >= size(res,2)
            break
        end
    end

    if isempty(high) == 0
        high_cell{1,h} = high;
        h = h + 1; 
    end
end

%% Plot the results to check:

% plot stochastic results
plot(res(1,:), res(2,:), 'Color',[0.4980,0.4980,0.4980])
hold on
% Plot low states:
low_tt = zeros(1,size(low_cell,2));
for i = 1:size(low_cell,2)
    time =  low_cell{1,i}(1,:);
    h = ones(size(time))*H_low_ss;
    if length(time) > 100
        low_tt(i) = time(end) - time(1);
        plot(time, h, 'Color',[0,0.45,0.74], 'LineWidth',2)
    end
end

% Plot high states:
high_tt = zeros(1,size(high_cell,2));
for i = 1:size(high_cell,2)
    time =  high_cell{1,i}(1,:);
    h = ones(size(time))*H_high_ss;
    if length(time) >100
        high_tt(i) = time(end) - time(1);
        plot(time, h, 'Color',[1,0,0], 'LineWidth',2)
    end
end

%% Save the figure
fig_title = "TypeI_mdl1_stochastic_NEW_long_r0"+ num2str(r*100) + "_S0" + num2str(S*100) + "_KH" + num2str(K_H) + "_" + sim_num;
fpath = "C:\Users\marko\OneDrive\Documents\Fakultet\06 Doktorske Studije\TA_typeI_models\stochastic_trajectories_figures";
saveas(gcf, fullfile(fpath, fig_title), 'jpeg')
saveas(gcf, fullfile(fpath, fig_title), "fig")
%% Calculate average transition time:

fprintf("Low toxin state duration:\nmin = %.2f min\nmax= %.2f min\nmean= %.2f min\n\n",min(low_tt),max(low_tt), mean(low_tt))

fprintf("High toxin state duration:\nmin = %.2f min \nmax= %.2f min \nmean= %.2f min \n\n",min(high_tt),max(high_tt), mean(high_tt))

tt_descriptives = table([min(low_tt);max(low_tt); mean(low_tt); median(low_tt)],[min(high_tt); max(high_tt); mean(high_tt); median(high_tt)],'VariableNames',["Low Toxin State tt [min]", "High Toxin State tt [min]"]);
tt_descriptives.Row = {'min','max','mean','median'};

%% Save the results:

filename = fig_title + '.mat';

save(filename, "tt_descriptives", "low_tt", "high_tt","res")

