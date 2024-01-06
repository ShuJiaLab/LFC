function histdata = plothist(data0,data30,data60,data120,data300)
%PLOTHIST Summary of this function goes here
%   Detailed explanation goes here
[dist0count,dist0edges] = histcounts(data0,30);
dist0edges = (dist0edges(1:end-1) + dist0edges(2:end))/2;
dist0count = dist0count/max(dist0count);

[dist30count,dist30edges] = histcounts(data30,30);
dist30edges = (dist30edges(1:end-1) + dist30edges(2:end))/2;
dist30count = dist30count/max(dist30count);

[dist60count,dist60edges] = histcounts(data60,30);
dist60edges = (dist60edges(1:end-1) + dist60edges(2:end))/2;
dist60count = dist60count/max(dist60count);

[dist120count,dist120edges] = histcounts(data120,20);
dist120edges = (dist120edges(1:end-1) + dist120edges(2:end))/2;
dist120count = dist120count/max(dist120count);

[dist300count,dist300edges] = histcounts(data300,20);
dist300edges = (dist300edges(1:end-1) + dist300edges(2:end))/2;
dist300count = dist300count/max(dist300count);

histdata = {[dist0edges',dist0count'];
            [dist30edges',dist30count'];
            [dist60edges',dist60count'];
            [dist120edges',dist120count'];
            [dist300edges',dist300count']};

figure(5),clf,
plot(dist0edges,dist0count,'red'),hold on;
plot(dist30edges,dist30count,'yellow'),hold on;
plot(dist60edges,dist60count,'blue'),hold on;
plot(dist120edges,dist120count,'green'),hold on;
plot(dist300edges,dist300count,'cyan'),hold on;
legend(["0 min","30 min","60 min","120 min","300 min"]);
end

