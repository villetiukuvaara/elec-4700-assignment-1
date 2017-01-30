%% Part 3: Enchancements

% The non-periodic top and bottom boundaries can be set to be either
% specular (1) or diffusive (0)
top_specular = 0;
bottom_specular = 0;

% Once again, the probability of scattering is
p_scat = 1 - exp(-time_step/0.2e-12);

% Each row corresponds to an electron with the position and velocities
% [x y vx vy]
state = zeros(population_size, 4);
trajectories = zeros(iterations, plot_population*2);
temperature = zeros(iterations,1);

% The vertices of the bottom-left and top-right corners of the boxes
% Also, each box can seperately be set to be specular or diffusive
% 1 = specular
% 0 = diffusive
boxes = 1e-9.*[80 120 0 40; 80 120 60 100];
boxes_specular = [0 1];

v_pdf = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/m));

% Generate an initial population
for i = 1:population_size
    angle = rand*2*pi;
    state(i,:) = [length*rand height*rand random(v_pdf) random(v_pdf)];
    
    % Make sure no particles start in a box
    while(in_box(state(i,1:2), boxes))
        state(i,1:2) = [length*rand height*rand];
    end
end

for i = 1:iterations
    %Update positions
    state(:,1:2) = state(:,1:2) + time_step.*state(:,3:4);
    
    j = state(:,1) > length;
    state(j,1) = state(j,1) - length;
    
    j = state(:,1) < 0;
    state(j,1) = state(j,1) + length;
    
    j = state(:,2) > height;

    if(top_specular)
        state(j,2) = 2*height - state(j,2);
        state(j,4) = -state(j,4);
    else % Re-thermaized
        % The electron bounces off at a random angle
        state(j,2) = height;
        v = sqrt(state(j,3).^2 + state(j,4).^2);
        angle = rand([sum(j),1])*2*pi;
        state(j,3) = v.*cos(angle);
        state(j,4) = -abs(v.*sin(angle));
    end
    
    j = state(:,2) < 0;
    
    if(bottom_specular)
        state(j,2) = -state(j,2);
        state(j,4) = -state(j,4);
    else % Re-thermalized
        % The electron bounces off at a random angle
        state(j,2) = 0;
        v = sqrt(state(j,3).^2 + state(j,4).^2);
        angle = rand([sum(j),1])*2*pi;
        state(j,3) = v.*cos(angle);
        state(j,4) = abs(v.*sin(angle));
    end
    
    % Look for particles that have "entered" a box and move them to
    % where they should be. Note that boxes should not overlap for this
    % to work
    for j=1:population_size
        box_num = in_box(state(j,1:2), boxes);
        while(box_num ~= 0)
            x_dist = 0;
            new_x = 0;
            if(state(j,3) > 0)
                x_dist = state(j,1) - boxes(box_num,1);
                new_x = boxes(box_num,1);
            else
                x_dist = boxes(box_num,2) - state(j,1);
                new_x = boxes(box_num,2);
            end
            
            y_dist = 0;
            new_y = 0;
            if(state(j,4) > 0)
                y_dist = state(j,2) - boxes(box_num, 3);
                new_y = boxes(box_num, 3);
            else
                y_dist = boxes(box_num, 4) - state(j,2);
                new_y = boxes(box_num, 4);
            end
            
            if(x_dist < y_dist)
                state(j,1) = new_x;
                if(~boxes_specular(box_num))
                    sgn = -sign(state(j,3));
                    v = sqrt(state(j,3).^2 + state(j,4).^2);
                    angle = rand()*2*pi;
                    state(j,3) = sgn.*abs(v.*cos(angle));
                    state(j,4) = v.*sin(angle);
                else
                    state(j,3) = -state(j,3);
                end
            else
                state(j,2) = new_y;
                if(~boxes_specular(box_num))
                    sgn = -sign(state(j,4));
                    v = sqrt(state(j,3).^2 + state(j,4).^2);
                    angle = rand()*2*pi;
                    state(j,3) = v.*cos(angle);
                    state(j,4) = sgn.*abs(v.*sin(angle));
                else
                    state(j,4) = -state(j,4);
                end
            end
            
            box_num = in_box(state(j,1:2), boxes);
        end
    end
    
    
    % Scatter particles
    j = rand(population_size, 1) < p_scat;
    state(j,3:4) = random(v_pdf, [sum(j),2]);
    
    % Record the temperature
    temperature(i) = (sum(state(:,3).^2) + sum(state(:,4).^2))*m/k/2/population_size;
    
    % Record positions for subset of particles that will be graphed
    % if show_movie = false
    for j=1:plot_population
        trajectories(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    
    % Update the movie every 5 iterations
    if show_movie && mod(i,5) == 0
        figure(3);
        subplot(3,1,1);
        hold off;
        plot(state(1:plot_population,1)./1e-9, state(1:plot_population,2)./1e-9, 'o');
        hold on;
        
        % Plot the boxes
        for j=1:size(boxes,1)
           plot([boxes(j, 1) boxes(j, 1) boxes(j, 2) boxes(j, 2) boxes(j, 1)]./1e-9,...
               [boxes(j, 3) boxes(j, 4) boxes(j, 4) boxes(j, 3) boxes(j, 3)]./1e-9, 'k-');
        end
        
        axis([0 length/1e-9 0 height/1e-9]);
        title(sprintf('Trajectories for %d of %d Electrons (Part 3)',...
        plot_population, population_size));
        xlabel('x (nm)');
        ylabel('y (nm)');
        if i > 1
            subplot(3,1,2);
            hold off;
            plot(time_step*(0:i-1), temperature(1:i));
            axis([0 time_step*iterations min(temperature(1:i))*0.98 max(temperature)*1.02]);
            title('Semiconductor Temperature');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
        end
        subplot(3,1,3);
        v = sqrt(state(:,3).^2 + state(:,4).^2);
        title('Histogram of Electron Speeds');
        histogram(v);
        xlabel('Speed (m/s)');
        ylabel('Number of particles');
        
        pause(0.05);
    end
end

% Show trajectories after the movie is over
figure(3);
subplot(3,1,1);
title(sprintf('Electron Trajectories for %d of %d Electrons (Part 3)',...
    plot_population, population_size));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 length/1e-9 0 height/1e-9]);
hold on;
for i=1:plot_population
    plot(trajectories(:,i*2)./1e-9, trajectories(:,i*2+1)./1e-9, '.');
    
end

% Plot the boxes
for j=1:size(boxes,1)
   plot([boxes(j, 1) boxes(j, 1) boxes(j, 2) boxes(j, 2) boxes(j, 1)]./1e-9,...
       [boxes(j, 3) boxes(j, 4) boxes(j, 4) boxes(j, 3) boxes(j, 3)]./1e-9, 'k-');
end

if(~show_movie)
    subplot(3,1,2);
    hold off;
    plot(time_step*(0:iterations-1), temperature);
    axis([0 time_step*iterations min(temperature)*0.98 max(temperature)*1.02]);
    title('Semiconductor Temperature');
    xlabel('Time (s)');
    ylabel('Temperature (K)');
end

subplot(3,1,3);
v = sqrt(state(:,3).^2 + state(:,4).^2);
title('Histogram of Electron Speeds');
histogram(v);
xlabel('Speed (m/s)');
ylabel('Number of particles');

% Create electron density map
density = hist3(state(:,1:2),[200 100])';

% Smooth out the electron density map
N = 20;
sigma = 3;
[x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(4);
subplot(1,2,1);
imagesc(conv2(density,f,'same'));
set(gca,'YDir','normal');
title('Electron Density');
xlabel('x (nm)');
ylabel('y (nm)');

% Create temperature map
temp_sum_x = zeros(ceil(length/1e-9),ceil(height/1e-9));
temp_sum_y = zeros(ceil(length/1e-9),ceil(height/1e-9));
temp_num = zeros(ceil(length/1e-9),ceil(height/1e-9));

for i=1:population_size
    x = floor(state(i,1)/1e-9);
    y = floor(state(i,2)/1e-9);
    if(x==0)
        x = 1;
    end
    if(y==0)
        y= 1;
    end
    temp_sum_y(x,y) = temp_sum_y(x,y) + state(i,3)^2;
    temp_sum_x(x,y) = temp_sum_x(x,y) + state(i,4)^2;
    temp_num(x,y) = temp_num(x,y) + 1;
end

% Create 2D histogram of temperatures over the region
temp = (temp_sum_x + temp_sum_y).*m./k./2./temp_num;
temp(isnan(temp)) = 0;
temp = temp';

% Smooth out the histogram (temperature map)
N = 20;
sigma = 3;
[x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
subplot(1,2,2);
imagesc(conv2(temp,f,'same'));
set(gca,'YDir','normal');
title('Temperature Map');
xlabel('x (nm)');
ylabel('y (nm)');