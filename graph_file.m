
function f = graph_file()

labRuntime = 'Running time (seconds)';
labExpansions = 'Number of expanded states';
labSolution = 'Shortest solution';



% EXAMPLE #1 OF HOW TO GRAPH STUFF
% Change the following if necessary
filegroup = 'lame';
xvals = [1, 2, 4, 8, 16];
xlab = 'Number of nodes/heuristics (one machine)';
titleRuntime = 'Running time with varying number of nodes';
titleExpansions = 'States expanded with varying number of nodes';
titleSolution = 'Solution length with varying number of nodes';

% No need to change anything from this point on
[runtime, expansions, solution, error_r, error_e, error_s] = get_data(filegroup, xvals);

fig_handle = graph_with_error_bars(xvals, runtime, error_r, xlab, labRuntime, true, true, titleRuntime);
print(fig_handle, ['graphs/', filegroup, '_runningtime.png'], '-dpng', '-r0');
fig_handle = graph_with_error_bars(xvals, expansions, error_e, xlab, labExpansions, true, true, titleExpansions);
print(fig_handle, ['graphs/', filegroup, '_expansions.png'], '-dpng', '-r0');
fig_handle = graph_with_error_bars(xvals, solution, error_s, xlab, labSolution, true, false, titleSolution);
print(fig_handle, ['graphs/', filegroup, '_solutionquality.png'], '-dpng', '-r0');
% /EXAMPLE #1



% EXAMPLE #2 OF HOW TO GRAPH STUFF
% Change the following if necessary
filegroup1 = 'lame';
filegroup2 = 'cool';
xvals = [1, 2, 4, 8, 16];
xlab = 'Number of nodes/heuristics (one machine)';
titleRuntime = 'Running time with varying number of nodes';
titleExpansions = 'States expanded with varying number of nodes';
titleSolution = 'Solution length with varying number of nodes';
legendlab1 = 'Lame SMHA*';
legendlab2 = 'Cool SMHA*';

% Do not need to change the following (I think)
[runtime1, expansions1, solution1, error_r1, error_e1, error_s1] = get_data(filegroup1, xvals);
[runtime2, expansions2, solution2, error_r2, error_e2, error_s2] = get_data(filegroup2, xvals);

fig_handle = multigraph_with_error_bars(xvals, runtime1, runtime2, error_r1, error_r2, xlab, labRuntime, legendlab1, legendlab2, true, true, titleRuntime);
print(fig_handle, ['graphs/', filegroup1, filegroup2, 'comp_runningtime.png'], '-dpng', '-r0');
fig_handle = multigraph_with_error_bars(xvals, expansions1, expansions2, error_e1, error_e2, xlab, labExpansions, legendlab1, legendlab2, true, true, titleExpansions);
print(fig_handle, ['graphs/', filegroup1, filegroup2, 'comp_expansions.png'], '-dpng', '-r0');
fig_handle = multigraph_with_error_bars(xvals, solution1, solution2, error_s1, error_s2, xlab, labSolution, legendlab1, legendlab2, true, false, titleSolution);
print(fig_handle, ['graphs/', filegroup1, filegroup2, 'comp_solutionquality.png'], '-dpng', '-r0');
% /EXAMPLE #2





% EXAMPLE #3 FREQUENCY
filegroup = 'freq';
xvals = [100, 200, 400, 800, 1600, 3200, 6400, 12800];
xlab = 'Number of iterations between sync-ups (4 nodes on 2 machines)';
titleRuntime = 'Running time with varying communication frequency';
titleExpansions = 'States expanded with varying communication frequency';
titleSolution = 'Solution length with varying communication frequency';

% No need to change anything from this point on
[runtime, expansions, solution, error_r, error_e, error_s] = get_data(filegroup, xvals);

fig_handle = graph_with_error_bars(xvals, runtime, error_r, xlab, labRuntime, true, true, titleRuntime);
print(fig_handle, ['graphs/', filegroup, '_runningtime.png'], '-dpng', '-r0');
fig_handle = graph_with_error_bars(xvals, expansions, error_e, xlab, labExpansions, true, true, titleExpansions);
print(fig_handle, ['graphs/', filegroup, '_expansions.png'], '-dpng', '-r0');
fig_handle = graph_with_error_bars(xvals, solution, error_s, xlab, labSolution, true, false, titleSolution);
print(fig_handle, ['graphs/', filegroup, '_solutionquality.png'], '-dpng', '-r0');
% /EXAMPLE #3


% MORE EXAMPLES
%xvals = [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800];
%xlab = 'Number of iterations between sync-ups (4 nodes on 2 machines)';
%xlab = 'Number of machines (1 core each, syncing every 1000 iterations)';



end


function [runtime, expansions, solution, error_r, error_e, error_s] = get_data(filegroup, filelist)

numfiles = size(filelist, 2);
runtime = zeros(1, numfiles);
expansions = zeros(1, numfiles);
solution = zeros(1, numfiles);
error_r = zeros(2, numfiles);
error_e = zeros(2, numfiles);
error_s = zeros(2, numfiles);

for i = 1:numfiles
  filename = ['stats_', filegroup, '_', num2str(filelist(i)), '.csv'];
  M = dlmread(filename, ' ');

  r_sort = sort(M(:, 3));
  e_sort = sort(M(:, 4));
  s_sort = sort(M(:, 5));

  runtime(i) = mean(r_sort(2:end-1)); % average values over trials and store
  error_r(:,i) = [r_sort(1); r_sort(end)]; % save lowest and highest for error bars
  expansions(i) = mean(e_sort(2:end-1));
  error_e(:,i) = [e_sort(1); e_sort(end)];
  solution(i) = mean(s_sort(2:end-1));
  error_s(:,i) = [s_sort(1); s_sort(end)];
end

end



function fig_handle = graph_with_error_bars(xvals, yvals, err, xlab, ylab, logx, logy, graphtitle)

fig_handle = figure;
set(gcf,'PaperPositionMode','auto');

% Plot values
hold on;
eb_handle = errorbar(xvals, yvals, err(1,:), err(2,:));
data_handle = plot(xvals, yvals);

xlabel(xlab);
ylabel(ylab);
title(graphtitle);

% Make figure bigger
set(fig_handle, 'Position', [0 0 600 300]);

% Make lines prettier
set(eb_handle, 'Color', [.9 .5 1]);
set(data_handle, 'Color', [.6 .2 .8], 'Marker', 's', 'MarkerSize', 8, ...
                 'LineWidth', 4, 'MarkerFaceColor', [.8 0 .6], 'MarkerEdgeColor', [.8 0 .6]);

if logx
  set(gca,'XScale','log');

  % Use nicer x-axis labels
  set(gca,'XTick', xvals);

  % Set nicer x-axis bounds
  xlim([xvals(1)/1.1, xvals(end)*1.1]);
end


if logy
  set(gca,'YScale','log')

  % Use nicer y-axis labels
  yticks = get(gca, 'YTick');
  ytickStr = [];
  for ytick = yticks
    if (ytick < 0.1)
      ytickStr = [ytickStr, sprintf('%1.3f|',ytick)];
    elseif (ytick < 1)
      ytickStr = [ytickStr, sprintf('%1.2f|',ytick)];
    else
      ytickStr = [ytickStr, sprintf('%1.0f|',ytick)];
    end
  end
  set(gca,'YTickLabel',ytickStr);
end
    
hold off;

end







function fig_handle = multigraph_with_error_bars(xvals, yvals1, yvals2, err1, err2, xlab, ylab, legendlab1, legendlab2, logx, logy, graphtitle)

fig_handle = figure;
set(gcf,'PaperPositionMode','auto');

% Plot values
hold on;
eb_handle1 = errorbar(xvals, yvals1, err1(1,:), err1(2,:));
eb_handle2 = errorbar(xvals, yvals2, err2(1,:), err2(2,:));
data_handle1 = plot(xvals, yvals1);
data_handle2 = plot(xvals, yvals2);

xlabel(xlab);
ylabel(ylab);
legend([data_handle1, data_handle2], {legendlab1, legendlab2});
title(graphtitle);

% Make figure bigger
set(fig_handle, 'Position', [0 0 600 300]);

% Make lines prettier
set(eb_handle1, 'Color', [1 .8 .6]);
set(data_handle1, 'Color', [1 .4 .3], 'Marker', 's', 'MarkerSize', 8, ...
                 'LineWidth', 4, 'MarkerFaceColor', [1 .4 .3], 'MarkerEdgeColor', [1 .4 .3]);
set(eb_handle2, 'Color', [.6 .8 1]);
set(data_handle2, 'Color', [0 .7 .9], 'Marker', 's', 'MarkerSize', 8, ...
                 'LineWidth', 4, 'MarkerFaceColor', [0 .7 .9], 'MarkerEdgeColor', [0 .7 .9]);

if logx
  set(gca,'XScale','log');

  % Use nicer x-axis labels
  set(gca,'XTick', xvals);

  % Set nicer x-axis bounds
  xlim([xvals(1)/1.1, xvals(end)*1.1]);
end


if logy
  set(gca,'YScale','log');

  % Use nicer y-axis labels
  yticks = get(gca, 'YTick');
  ytickStr = [];
  for ytick = yticks
    if (ytick < 0.1)
      ytickStr = [ytickStr, sprintf('%1.3f|',ytick)];
    elseif (ytick < 1)
      ytickStr = [ytickStr, sprintf('%1.2f|',ytick)];
    else
      ytickStr = [ytickStr, sprintf('%1.0f|',ytick)];
    end
  end
  set(gca,'YTickLabel',ytickStr);
end
    
hold off;

end


