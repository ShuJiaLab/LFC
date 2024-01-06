% clear,clc
% dist0 = [];
% dist30 = [];
% dist60 = [];
% dist120 = [];
% dist300 = [];
%%
[dist0count,dist0edges] = histcounts(dists0,30);
dist0edges = (dist0edges(1:end-1) + dist0edges(2:end))/2;
% dist0count = dist0count/length(dists0);
dist0count = dist0count/max(dist0count);

[dist30count,dist30edges] = histcounts(dists30,30);
dist30edges = (dist30edges(1:end-1) + dist30edges(2:end))/2;
% dist30count = dist30count/length(dists30);
dist30count = dist30count/max(dist30count);

[dist60count,dist60edges] = histcounts(dists60,30);
dist60edges = (dist60edges(1:end-1) + dist60edges(2:end))/2;
% dist60count = dist60count/length(dists60);
dist60count = dist60count/max(dist60count);

[dist120count,dist120edges] = histcounts(dists120,20);
dist120edges = (dist120edges(1:end-1) + dist120edges(2:end))/2;
% dist120count = dist120count/length(dists120);
dist120count = dist120count/max(dist120count);

[dist300count,dist300edges] = histcounts(dists300,20);
dist300edges = (dist300edges(1:end-1) + dist300edges(2:end))/2;
% dist300count = dist300count/length(dists300);
dist300count = dist300count/max(dist300count);

dists = {dist0edges',dist0count',dist30edges',dist30count',dist60edges',dist60count',dist120edges',dist120count',dist300edges',dist300count'};
%%
f0 = fit(dist0edges',dist0count','gauss1');
f30 = fit(dist30edges',dist30count','gauss1');
f60 = fit(dist60edges',dist60count','gauss1');
f120 = fit(dist120edges',dist120count','gauss1');
f300 = fit(dist300edges',dist300count','gauss1');
%%
clc
figure(5),clf,
plot(dist0edges,dist0count,'red'),hold on;
plot(dist30edges,dist30count,'yellow'),hold on;
plot(dist60edges,dist60count,'blue'),hold on;
plot(dist120edges,dist120count,'green'),hold on;
plot(dist300edges,dist300count,'cyan'),hold on;
legend(["0 min","30 min","60 min","120 min","300 min"]);
%%
clc
figure(5),
plot(f0,'red'),hold on;
plot(f30,'yellow'),hold on;
plot(f60,'blue'),hold on;
plot(f120,'green'),hold on;
plot(f300,'cyan'),hold on;
legend(["0 min"," ","30 min"," ","60 min"," ","120 min"," ","300 min"," "]);