


function output_ax = dual_plot_ini(input, var_1, var_2, do_log_1, do_log_2, i, cmap)

output_ax = struct();

tmp_var_1 = var_1( (input.XX < input.x_range(2)) & (input.XX > input.x_range(1)) );
tmp_var_2 = var_2( (input.XX < input.x_range(2)) & (input.XX > input.x_range(1)) );

figure(input.fig);
output_ax.ax_11 = subplot(211);
set(gca, 'fontsize', input.fontsize);
colormap(input.(char(cmap)));
if do_log_1
    var_1(var_1<1e-30)=1e-30;
    tmp_var_1 = log10(var_1( (input.XX < input.x_range(2)) & (input.XX > input.x_range(1)) ));
    c_max = max(tmp_var_1(:));
    tmp_var_1(tmp_var_1==-30) = c_max;
    c_min = min(tmp_var_1(:));
    if c_min==c_max; c_min = c_max - 1; c_max = c_max + 1; end;
    pcolor(input.ZZ,input.XX,log10(var_1)), shading flat, cb = colorbar(); set(cb, 'fontsize', input.fontsize);
    caxis([c_min c_max]);
    xlabel(' z (um) '), ylabel(' x (um) '), output_ax.title_1 = title(input.title{1}{2});
else
    pcolor(input.ZZ,input.XX,var_1), shading flat, cb = colorbar(); set(cb, 'fontsize', input.fontsize);
    c_min = min(tmp_var_1(:));
    c_max = max(tmp_var_1(:));
    if c_min==c_max; c_min = c_max - 1; c_max = c_max + 1; end;
    if strcmp(cmap, 'bwr'); c_max = max([abs(c_min) abs(c_max)]); c_min = -c_max; end;
    caxis([c_min c_max]);
    xlabel(' z (um) '), ylabel(' x (um) '), output_ax.title_1 = title(input.title{1}{1});
end
ylim(input.x_range);
output_ax.ax_12 = get(output_ax.ax_11, 'Children');

output_ax.ax_21 = subplot(212);
set(gca, 'fontsize', input.fontsize);
colormap(input.(char(cmap)));
if do_log_2
    var_2(var_2<1e-30)=1e-30;
    tmp_var_2 = log10(var_2( (input.XX < input.x_range(2)) & (input.XX > input.x_range(1)) ));
    c_max = max(tmp_var_2(:));
    tmp_var_2(tmp_var_2==-30) = c_max;
    c_min = min(tmp_var_2(:));
    if c_min==c_max; c_min = c_max - 1; c_max = c_max + 1; end;
    pcolor(input.ZZ,input.XX,log10(var_2)), shading flat, cb = colorbar(); set(cb, 'fontsize', input.fontsize);
    caxis([c_min c_max]);
    xlabel(' z (um) '), ylabel(' x (um) '), output_ax.title_2 = title(input.title{2}{2});
else
    pcolor(input.ZZ,input.XX,var_2), shading flat, cb = colorbar(); set(cb, 'fontsize', input.fontsize);
    c_min = min(tmp_var_2(:));
    c_max = max(tmp_var_2(:));
    if c_min==c_max; c_min = c_max - 1; c_max = c_max + 1; end;
    if strcmp(cmap, 'bwr'); c_max = max([abs(c_min) abs(c_max)]); c_min = -c_max; end;
    caxis([c_min c_max]);
    xlabel(' z (um) '), ylabel(' x (um) '), output_ax.title_2 = title(input.title{2}{1});
end
ylim(input.x_range);
output_ax.ax_22 = get(gca, 'Children');

if input.do_save
    filename = [input.path 'frame_' num2str(i, '%.5d') '.png'];
    saveas(input.fig, filename, 'png');
else
    pause(0.001);
end


    
end