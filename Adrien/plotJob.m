function plotJob(njob,nmachine,type)
close all;

% read input lines
prefix=strcat('results',num2str(njob),'x',num2str(nmachine),'_',type);
f=fopen(strcat('./output/',prefix,'.out'));
tline = fgetl(f);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(f);
end
fclose(f);


n = sscanf(tlines{1}, '%d %d'); %number of jobs and maintenance tasks
OF =  sscanf(tlines{2}, '%f %f %f %f'); %Factors and function objective value
OFVS =  sscanf(tlines{3}, '%f %f %f'); %Objective values
data_jobs = zeros(n(1), 5); %scheduling data of jobs
data_maint = zeros(n(2), 4); %scheduling data of maintenances tasks
for i=1:n(1)
  data_jobs(i,:) = sscanf(tlines{i+3}, '%d %d %d %f %f');
end
for i=1:n(2)
  data_maint(i,:) = sscanf(tlines{i+n(1)+3}, '%d %d %d %f %f');
end

colors = lines(n(1));
figure;
%Plot the bloc of each operation
for i=1:n(1)
  X = data_jobs(i,4);
  Y = data_jobs(i,3);
  W = data_jobs(i,5) - data_jobs(i,4);
  H = 1;
  color = colors(data_jobs(i,1)+1,:) + 0.1*floor((data_jobs(i,1)+1)/8)*ones(1,3);
  color = mod(color,1);
  rectangle('Position',[X Y W H],'FaceColor',color);
  text(X,Y+H/2, ['O_{',num2str(data_jobs(i,1)+1),',',num2str(data_jobs(i,2)+1),'}']);
end
%Plot the bloc of each maintenance task
for i=1:n(2)
  X = data_maint(i,3);
  Y = data_maint(i,1);
  W = data_maint(i,4) - data_maint(i,3);
  H = 1;
  color = [0.7 0.7 0.7];
  rectangle('Position',[X Y W H],'FaceColor',color);
  text(X,Y+H/2, ['PM_{',num2str(data_maint(i,1)+1),',',num2str(data_maint(i,2)+1),'}']);
end
maintitle = ['W_{tot}=',num2str(OFVS(1)),'  W_{max}=',num2str(OFVS(2)),'  C_{max}=',num2str(OFVS(3))];
subtitle = ['F=',num2str(OF(1)),'*W_{tot}+',num2str(OF(2)),'*W_{max}+',num2str(OF(3)),'*C_{max}','=',num2str(OF(4))];
title ( [maintitle,' ',subtitle]);
saveas(gcf,strcat('./img/',prefix,'.png'));
