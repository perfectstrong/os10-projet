vals=[4 5;8 8;10 10;15 10];
close all;
for i=1:length(vals)
  plotFront(vals(i,1),vals(i,2));
end
