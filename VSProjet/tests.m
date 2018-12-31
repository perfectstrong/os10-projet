types={'Cmax','Wmax','Wtot','F050203','F050302'};
vals=[4 5;8 8;10 10;15 10];
for i=1:length(vals)
    for j=1:length(types)
        plotJob(vals(i,1),vals(i,2),types{j});
    end
end
