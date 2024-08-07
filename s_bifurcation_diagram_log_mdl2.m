function zero_vec = s_bifurcation_diagram_log_mdl2(s_vec,r, axis_font, label_font, colormap,lw)
%r_bifurcation_diagram draws H vs r_vec bifurcation diagram
%   This function draws a bifurcation diagram for the values
%   of parameter r given in r_vec and a set value of parameter s, with H on
%   log scale:

zero_vec = [];
for s = s_vec
    pol = [ 1 , 1-1/s , 0, 1/(r*s) , 1 + 1/(r*s) , 1-1/s ];
    zerov = roots(pol);
    zero_vec = [zero_vec; zerov'];
end

stability_vec = [];
lena = size(zero_vec,1);
s_high_vec = [];
high_zero_vec = [];
s_unstable_vec = [];
unstable_zero_vec = [];
s_low_vec = [];
low_zero_vec = [];

for i = 1:lena
    s = s_vec(1,i);
    zerov = zero_vec( i , : );
    root_1 = zerov(1);
    root_2 = zerov(2);
    root_3 = zerov(3);
    root_4 = zerov(4);
    root_5 = zerov(5);

    if isreal( root_1 )
        high_zero_vec = [ high_zero_vec  root_1 ];
        s_high_vec = [ s_high_vec s ];
    end


    if isreal( root_2 )
        unstable_zero_vec = [ unstable_zero_vec  root_2 ];
        s_unstable_vec = [ s_unstable_vec s ];
    end


    if isreal( root_4 )
        unstable_zero_vec = [ unstable_zero_vec  root_4 ];
        s_unstable_vec = [  s_unstable_vec   s ];
    end

    if isreal( root_5 )
        low_zero_vec = [ low_zero_vec  root_5 ];
        s_low_vec = [  s_low_vec   s ];
    end

    sum = isreal( root_1 ) + isreal( root_2 ) + isreal( root_3 )  + isreal( root_4 ) + isreal( root_5 ) ;
    stability_vec = [ stability_vec  [ s ; sum ] ]   ;

end

plot(s_low_vec, log(low_zero_vec), Color=colormap(1,:), LineWidth=lw, LineStyle= "-")
hold on
plot(s_unstable_vec,log(unstable_zero_vec), Color=colormap(2,:), LineWidth=lw, LineStyle="--")
plot(s_high_vec, log(high_zero_vec), Color=colormap(3,:), LineWidth=lw, LineStyle="-")
ax = gca;
ax.FontSize = axis_font;
xlabel('S', 'FontSize',label_font)
hold off
end