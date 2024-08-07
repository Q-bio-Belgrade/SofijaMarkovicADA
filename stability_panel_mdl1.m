%% Drawing Model 1 combined stability diagram panel
% This code draws panel including both S vs r double bifurcation
% diagram (stability diagram) and six single bifurcation diagrams for S and r values indicated
% by lines in stability diagram. 

% Sofija Markovic
% 13/4/2024

clear

%% Set Figure Properties

% Stability diagram formatting:
line_label_fs = 18;
axis_fs = 22;
ax_label_fs = 32;
stability_S_line_col  = [112,48,100]/255;
stability_r_line_col = [242, 170, 132]/255; 
stability_main = [64, 64, 64]/255; 
stability_fill = [127, 127, 127]/255; 
r_s_lines_col = [89, 89, 89]/255;
stability_main_lw = 2;
stability_x_y_lw = 1;

% bifurcation diagrams formatting:
axis_font = 16;
label_font = 20;
annot_font = 16;
colormap = [0,0.45,0.74;0.47,0.67,0.19;1,0,0;]; 
lw = 1.5;

% Include titles to bif_diags:
bd_title = false;
bd_annot = true;
bd_label = true;

% Type of panel:
% 1: 4x4
% 2: 8x7
% 3: 7x8

panel = 2;

switch panel
    case 1
        tiledlayout(4,4, "TileSpacing","tight")
    case 2
        tiledlayout(8,7,"TileSpacing","tight")
    case 3
        tiledlayout(7,8, "TileSpacing","tight")
    otherwise
        error('Please state panel as integer 1, 2 or 3')
end

%% Calculations:

% set the values of H:
H = (1/3)^(0.25):0.0001:15;

% calculate s and r: 
s = (3*H.^4 - 1)./(4*H.^5 + 3*H.^4 - 1);
r = (4*H.^5 + 3*H.^4 - 1)./((1 + H.^4).^2);

% Define S and r vec for bif diags
s_vec =  0.0001:0.0001:0.8;
r_vec = 0.0001:0.0001:2;
% For fill function:
X = [0, r]; 
Y = [0, s];

%% S vs r stability diagram

switch panel
    case 1
        nexttile(1,[3,3])
    case 2
        nexttile(1,[6,5])
    case 3
        nexttile(1,[5,6])
end

fill(X,Y,stability_fill, 'EdgeColor','none','FaceAlpha',0.6)
hold on
plot(r,s, Color=stability_main, LineWidth=stability_main_lw)
ax1 = gca;
ax1.FontSize = axis_fs;
ax1.TickLength = [0.06 0.105]; % ovo je 3x default vrednost

% s values for r bifurcation diagram:
yl1 = yline(0.35,'--', 'S = 0.35',Color=r_s_lines_col,LineWidth=stability_x_y_lw);
yl1.FontSize = line_label_fs;
yl2 = yline(0.2,'-', 'S = 0.2', Color=r_s_lines_col,LineWidth=stability_x_y_lw);
yl2.FontSize = line_label_fs;
yl3 = yline(0.05,':', 'S = 0.05', Color=r_s_lines_col,LineWidth=stability_x_y_lw);
yl3.FontSize = line_label_fs;
% r values for s bifurcation diagram:
xl1 = xline(1.7,'--', 'r = 1.7', Color=r_s_lines_col,LineWidth=stability_x_y_lw);
xl1.FontSize = line_label_fs;
xl2 = xline(0.8,'-', 'r = 0.8', Color=r_s_lines_col,LineWidth=stability_x_y_lw);
xl2.FontSize = line_label_fs;
xl3 = xline(0.2,':', 'r = 0.2', Color=r_s_lines_col,LineWidth=stability_x_y_lw);
xl3.FontSize = line_label_fs;

xlabel('r', 'FontSize',ax_label_fs)
ylabel('S','FontSize', ax_label_fs)
box("on")
xlim([0,2])
ylim([0,0.4])
yticks([0, 0.1, 0.2, 0.3, 0.4])
xticks([0, 0.5, 1, 1.5, 2])

if bd_label == true
        range_x = max(xlim) - min(xlim);
        range_y = max(ylim) - min(ylim);
        pos_y = max(ylim) - range_y*0.1;
        pos_x = min(xlim) + range_x*0.02;
        text(pos_x,pos_y,'A', 'Horiz','left', 'Vert','bottom','FontSize', 25)
end
hold off

%% r diagram for S  = 0.35
switch panel 
    case 1
        nexttile
    otherwise
        nexttile([2,2])
end

r_bifurcation_diagram_mdl1(r_vec, 0.35, axis_font, label_font, colormap,lw, {'-','--','-'})
set(gca,'TickLength',[0.06 0.105]);

if bd_title == true
    title('S = 0.35')
end
if bd_annot == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = min(ylim) + range_y*0.05;
    pos_x = max(xlim) - range_x*0.05;
    text(pos_x, pos_y,'S = 0.35', 'Horiz','right', 'Vert','bottom','FontSize',annot_font)
end

if bd_label == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.05;
    pos_x = min(xlim) + range_x*0.075;
    text(pos_x, pos_y,'E', 'Horiz','right', 'Vert','bottom', 'FontSize',annot_font)
end
%% r diagram for S  = 0.2

switch panel 
    case 1
        nexttile
    otherwise
        nexttile([2,2])
end

r_bifurcation_diagram_mdl1(r_vec, 0.2, axis_font, label_font, colormap,lw,{'-','--','-'});
set(gca,'TickLength',[0.06 0.105])

if bd_title == true
    title('S = 0.2')
end
if bd_annot == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = min(ylim) + range_y*0.05;
    pos_x = max(xlim) - range_x*0.05;
    text(pos_x, pos_y,'S = 0.2', 'Horiz','right', 'Vert','bottom','FontSize',annot_font)
end

if bd_label == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.175;
    pos_x = min(xlim) + range_x*0.075;
    text(pos_x, pos_y,'F', 'Horiz','right', 'Vert','bottom', 'FontSize',annot_font)
end

ArrowDown = annotation('arrow','LineStyle','--','LineWidth',2);
ArrowDown.Parent = gca;
ArrowDown.X = [0.19,0.19]; 
ArrowDown.Y = [2.7 ,0.3];

ArrowUp = annotation('arrow','LineStyle','--','LineWidth',2);
ArrowUp.Parent = gca;
ArrowUp.X = [0.87,0.87];
ArrowUp.Y = [1.08 ,3.8];
%% r diagram for S  = 0.05
switch panel 
    case 1
        nexttile
    otherwise
        nexttile([2,2])
end

r_bifurcation_diagram_log_mdl1(r_vec, 0.05, axis_font, label_font, colormap,lw,{'-','--','-'});
set(gca,'TickLength',[0.06 0.105])
xlim([0,2])
ylim([-8,6])
yticks([-8, -1, 6])
xticks([0, 1, 2])
if bd_title == true
    title('S = 0.05')
end

if bd_annot == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = min(ylim) + range_y*0.05;
    pos_x = max(xlim) - range_x*0.05;
    text(pos_x, pos_y,'S = 0.05', 'Horiz','right', 'Vert','bottom','FontSize',annot_font)
end

if bd_label == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.175;
    pos_x = min(xlim) + range_x*0.075;
    text(pos_x, pos_y,'G', 'Horiz','right', 'Vert','bottom', 'FontSize',annot_font)
end

%% S diagram for r  = 0.2
switch panel 
    case 1
        nexttile
    otherwise
        nexttile([2,2])
end

s_bifurcation_diagram_log_mdl1(s_vec, 0.2, axis_font, label_font, colormap,lw);
set(gca,'TickLength',[0.06 0.105])
ylabel('log(H)', 'FontSize',label_font)

if bd_title == true
    title('r = 0.2')
end
if bd_annot == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.05;
    pos_x = max(xlim) - range_x*0.05;
    text(pos_x, pos_y,'r = 0.2', 'Horiz','right', 'Vert','top','FontSize',annot_font)
end

if bd_label == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.175;
    pos_x = min(xlim) + range_x*0.075;
    text(pos_x, pos_y,'B', 'Horiz','right', 'Vert','bottom', 'FontSize',annot_font)
end

%% S diagram for r  = 0.8
switch panel 
    case 1
        nexttile
    otherwise
        nexttile([2,2])
end

s_bifurcation_diagram_log_mdl1(s_vec, 0.8, axis_font, label_font, colormap,lw);
set(gca,'TickLength',[0.06 0.105])

if bd_title == true    
    title('r = 0.8')
end
if bd_annot == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.05;
    pos_x = max(xlim) - range_x*0.05;
    text(pos_x, pos_y,'r = 0.8', 'Horiz','right', 'Vert','top','FontSize',annot_font)
end

if bd_label == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.175;
    pos_x = min(xlim) + range_x*0.075;
    text(pos_x, pos_y,'C', 'Horiz','right', 'Vert','bottom', 'FontSize',annot_font)
end

%% S diagram for r  = 1.7
switch panel 
    case 1
        nexttile
    otherwise
        nexttile([2,2])
end
s_bifurcation_diagram_log_large_r_mdl1(s_vec, 1.7, axis_font, label_font, colormap,lw);
set(gca,'TickLength',[0.06 0.105])

if bd_title == true
    title('r = 1.7')
end
if bd_annot == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.05;
    pos_x = max(xlim) - range_x*0.05;
    text(pos_x, pos_y,'r = 1.7', 'Horiz','right', 'Vert','top','FontSize',annot_font)
end

if bd_label == true
    range_x = max(xlim) - min(xlim);
    range_y = max(ylim) - min(ylim);
    pos_y = max(ylim) - range_y*0.175;
    pos_x = min(xlim) + range_x*0.075;
    text(pos_x, pos_y,'D', 'Horiz','right', 'Vert','bottom', 'FontSize',annot_font)
end