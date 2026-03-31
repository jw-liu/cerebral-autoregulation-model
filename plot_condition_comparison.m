function plot_condition_comparison(results, net, params, c, fname)

mu_exp  = params.exp_mu(c, :); sig_exp = params.exp_sig(c, :); sim_cam = results.cam_rel(c, :);
if c == 1
    sim_nc = NaN(1, 10);
else
    sim_nc = results.nocam_rel(c, :);
end
idx_s = net.idx_source; idx_k = net.idx_sink;

if c == 2
    inlet_labels = {'ipsi-ICA', 'cont-ICA', 'ipsi-VA', 'cont-VA'};
    outlet_labels = {'ipsi-ACAA1', 'cont-ACA', 'ipsi-MCA', 'cont-MCA', 'ipsi-PCA', 'cont-PCA'};
elseif c == 3
    inlet_labels = {'ipsi-ICA', 'cont-ICA', 'ipsi-VA', 'cont-VA'};
    outlet_labels = {'hypoplastic ACAA1', 'cont-ACA', 'ipsi-MCA', 'cont-MCA', 'ipsi-PCA', 'cont-PCA'};
else
    inlet_labels = net.labels_10(idx_s);
    outlet_labels = net.labels_10(idx_k);
end

mu_plot  = [mu_exp(idx_s), NaN, mu_exp(idx_k)];
sig_plot = [sig_exp(idx_s),NaN, sig_exp(idx_k)];
nc_plot  = [sim_nc(idx_s),NaN, sim_nc(idx_k)];
cam_plot = [sim_cam(idx_s),NaN, sim_cam(idx_k)];
lab_plot = [inlet_labels, {''}, outlet_labels];
n = length(mu_plot); has_nc = any(~isnan(nc_plot));
clr_cam= [0.8 0 0]; clr_nocam = [0.2 0.7 0.2];

f = figure('Position', [50 400 820 390], 'Color', 'w','Visible', 'off');
if c == 1
    ax = axes('Position', [0.10 0.18 0.85 0.72]);
elseif c == 2
   ax = axes('Position', [0.10 0.2 0.85 0.71]);
elseif c == 3
       ax = axes('Position', [0.10 0.27 0.85 0.66]);
end
bar(ax, 1:n, mu_plot, 'FaceColor', [0.6 0.6 0.6],'BarWidth', 0.65, 'EdgeColor', 'k', 'LineWidth', 1.2);
hold on;
errorbar(1:n, mu_plot, sig_plot, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 8);
for i = 1:n
    if ~isnan(cam_plot(i))
        plot(i, cam_plot(i), 'o', 'MarkerSize', 10, 'MarkerFaceColor', clr_cam, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    end
    if has_nc && ~isnan(nc_plot(i))
        plot(i, nc_plot(i), '^', 'MarkerSize', 10,'MarkerFaceColor', clr_nocam, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    end
end

hl = xline(length(idx_s) + 1, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
set(get(get(hl, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
xticks(1:n); xticklabels(lab_plot); xtickangle(30);
ylabel('Relative Flow rate (%)', 'FontSize', 13, 'FontWeight', 'bold');
if c == 2
    ttl = 'Fetal-type P1 segment of PCA';
elseif c == 3
    ttl = 'Hypoplastic A1 segment of ACA';
else
    ttl = params.cond_titles{c};
end
title(ttl, 'FontWeight', 'bold', 'FontSize', 14); grid on; box on;
set(gca, 'FontSize', 13, 'FontName', 'Arial', 'LineWidth', 1.2, 'TickDir', 'out');
%ylim([0, max(mu_plot + sig_plot, [], 'omitnan') + 8]);
xlim([0 12]); ylim([0,45]);
if c == 3
    ylim([0,63]);
end

hB = bar(NaN, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k');
hE = errorbar(-99, 0, 1, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 8);
if has_nc
    hN = plot(NaN, NaN, '^', 'MarkerSize', 9, 'MarkerFaceColor', clr_nocam, 'MarkerEdgeColor', 'k');
    hC = plot(NaN, NaN, 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', clr_cam, 'MarkerEdgeColor', 'k');
    legend([hB, hE, hC, hN], {'Clinical mean', '\pm Standard Deviation','Simulation (with autoregulation)',...
        'Simulation (without autoregulation)'}, 'NumColumns', 1, 'Location', 'best', 'FontSize', 10);
else
    hC = plot(NaN, NaN, 'o', 'MarkerSize', 10,'MarkerFaceColor', clr_cam, 'MarkerEdgeColor', 'k');
    legend([hB, hE, hC], {'Clinical mean', '\pm Standard Deviation', 'Simulation (with autoregulation)'}, ...
        'Location', 'best', 'FontSize', 11);
end
saveas(f, [fname '.png']);
fprintf('  Saved: %s.png\n', fname);
end
