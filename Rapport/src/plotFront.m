function plotFront(njob,nmachine)
  % read input lines
  prefix=strcat('results',num2str(njob),'x',num2str(nmachine),'_ParetoFront');
  f=fopen(strcat('./output/',prefix,'.out'));
  tline = fgetl(f);
  tlines = cell(0,1);
  nb_points =0;
  while ischar(tline)
      tlines{end+1,1} = tline;
      tline = fgetl(f);
      nb_points++;
  end
  fclose(f);
  figure;
  colors = lines(16); % no more than 16 fronts
  data = zeros(nb_points,4);
  for i=1:nb_points
    data(i,:) = sscanf(tlines{i}, '%f %f %f %f');
  endfor
  %plot objective vector
  scatter3 (data(:,2),data(:,3),data(:,4), 50, colors(data(:,1)),'filled')
  text(data(:,2)+0.1,data(:,3)+0.1,data(:,4)+0.1,num2str(data(:,1)));
  %labels and title
  xlabel("W_{tot}");
  ylabel("W_{max}");
  zlabel("C_{max}");
  title(["Frontières de Pareto pour le problem ",num2str(njob)," x ", num2str(nmachine)]);
endfunction
