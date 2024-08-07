function zero_vec = r_bifurcation_diagram_log_mdl1(r_vec,s,axis_font, label_font, colormap,lw, ls)
% r_bifurcation_diagram draws H vs r_vec bifurcation diagram for type I TA
% model 1
%   This function draws a bifurcation diagram for the given value
%   of parameter s and a set of values of parameter r given in r_vec, with H on log scale:

zero_vec = [];
for r = r_vec
    pol = [ 1 , 1-1/s , 0, 0 , 1 + 1/(r*s) , 1-1/s ];
    zerov = roots(pol);
    zero_vec = [zero_vec; zerov'];
end

stability_vec = [];
lena = size(zero_vec,1);
r_high_vec = [];
high_zero_vec = [];
r_unstable_vec = [];
unstable_zero_vec = [];
r_low_vec = [];
low_zero_vec = [];

for i = 1:lena
    r = r_vec(1,i);
    zerov = zero_vec( i , : );
    root_1 = zerov(1);
    root_2 = zerov(2);
    root_3 = zerov(3);
    root_4 = zerov(4);
    root_5 = zerov(5);

    if isreal( root_1 )
        high_zero_vec = [ high_zero_vec  root_1 ];
        r_high_vec = [ r_high_vec r ];
    end


    if isreal( root_2 )
        unstable_zero_vec = [ unstable_zero_vec  root_2 ];
        r_unstable_vec = [ r_unstable_vec r ];
    end


    if isreal( root_4 )
        unstable_zero_vec = [ unstable_zero_vec  root_4 ];
        r_unstable_vec = [  r_unstable_vec   r ];
    end

    if isreal( root_5 )
        low_zero_vec = [ low_zero_vec  root_5 ];
        r_low_vec = [  r_low_vec   r ];
    end

    sum = isreal( root_1 ) + isreal( root_2 ) + isreal( root_3 )  + isreal( root_4 ) + isreal( root_5 ) ;
    stability_vec = [ stability_vec  [ r ; sum ] ]   ;

end

plot(r_low_vec, log(low_zero_vec), Color=colormap(1,:), LineWidth=lw, LineStyle= ls{1})
hold on
plot(r_unstable_vec,log(unstable_zero_vec), Color=colormap(2,:), LineWidth=lw, LineStyle=ls{2})
plot(r_high_vec, log(high_zero_vec), Color=colormap(3,:), LineWidth=lw, LineStyle=ls{3})
ax = gca;
ax.FontSize = axis_font;
xlabel('r', 'FontSize',label_font)
ylabel('log(H)', 'FontSize',label_font)

end