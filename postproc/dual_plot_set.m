


function dual_plot_set(input_ax, input, var_1, var_2, do_log_1, do_log_2, i, cmap)

tmp_var_1 = var_1( (input.XX < input.x_range(2)) & (input.XX > input.x_range(1)) );
tmp_var_2 = var_2( (input.XX < input.x_range(2)) & (input.XX > input.x_range(1)) );

if do_log_1
    var_1(var_1<1e-30)=1e-30;
    tmp_var_1 = log10(var_1( (input.XX < input.x_range(2)) & (input.XX > input.x_range(1)) ));
    c_max = max(tmp_var_1(:));
    tmp_var_1(tmp_var_1==-30) = c_max;
    c_min = min(tmp_var_1(:));
    if c_min==c_max; c_min = c_max - 1; c_max = c_max + 1; end;
    set(input_ax.ax_11, 'CLim', [c_min c_max]);
    set(input_ax.ax_12, 'CData', log10(var_1));
%     set(input_ax.ax_12(8), 'CData', log10(var_1));
    set(input_ax.title_1, 'String', input.title{1}{2});
else
    c_min = min(tmp_var_1(:));
    c_max = max(tmp_var_1(:));
    if c_min==c_max; c_min = c_max - 1; c_max = c_max + 1; end;
    if strcmp(cmap, 'bwr'); c_max = max([abs(c_min) abs(c_max)]); c_min = -c_max; end;
    set(input_ax.ax_11, 'CLim', [c_min c_max]);
    set(input_ax.ax_12, 'CData', var_1);
%     set(input_ax.ax_12(8), 'CData', var_1);
    set(input_ax.title_1, 'String', input.title{1}{1});
end
if do_log_2
    var_2(var_2<1e-30)=1e-30;
    tmp_var_2 = log10(var_2( (input.XX < input.x_range(2)) & (input.XX > input.x_range(1)) ));
    c_max = max(tmp_var_2(:));
    tmp_var_2(tmp_var_2==-30) = c_max;
    c_min = min(tmp_var_2(:));
    if c_min==c_max; c_min = c_max - 1; c_max = c_max + 1; end;
    set(input_ax.ax_21, 'CLim', [c_min c_max]);
    set(input_ax.ax_22, 'CData', log10(var_2));
%     set(input_ax.ax_22(8), 'CData', log10(var_2));
    set(input_ax.title_2, 'String', input.title{2}{2});
else
    c_min = min(tmp_var_2(:));
    c_max = max(tmp_var_2(:));
    if c_min==c_max; c_min = c_max - 1; c_max = c_max + 1; end;
    if strcmp(cmap, 'bwr'); c_max = max([abs(c_min) abs(c_max)]); c_min = -c_max; end;
    set(input_ax.ax_21, 'CLim', [c_min c_max]);
    set(input_ax.ax_22, 'CData', var_2);
%     set(input_ax.ax_22(8), 'CData', var_2);
    set(input_ax.title_2, 'String', input.title{2}{1});
end

if input.do_save
    filename = [input.path 'frame_' num2str(i, '%.5d') '.png'];
    saveas(input.fig, filename, 'png');
else
    pause(0.001);
end


end