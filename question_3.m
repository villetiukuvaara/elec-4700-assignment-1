%%Question 3

m0 = 9.10938356e-31;
m = 0.26*m0;
T = 300;
k = 1.38064852e-23;

%% Question 1
% To calculate the thermal energy, it is noted that there are two degrees
% of freedom for the electrons. Using Maxwell's principle of equipartition
% of energy,
%
% $$\overline{KE}=\frac{1}{2}kT = 2(\frac{1}{2}m\overline{v^2}) \Rightarrow \overline{v^2}=\frac{2kT}{m}$$

vth = sqrt(2*k*T/m);

%%
% Or 18.7 um/s. The mean free path, $l$, is simply

l = vth*0.2e-12;

%%
% Thus the mean free path is about 37.4 nm.

% Set up the simulation

width = 100e-9;
length = 200e-9;
population_size = 5000;
plot_population = 20;
time_step = width/vth/100;
iterations = 100;
show_movie = 1;
p_scat = 1 - exp(-time_step/0.2e-12);

% Each row corresponds to an electron with the position and velocities
% [x y vx vy]
state = zeros(population_size, 4);
trajectories = zeros(iterations, plot_population*2);
temperature = zeros(iterations,1);

% The vertices of the bottom-left and top-right corners of the boxes
boxes = 1e-9.*[80 120 0 40; 80 120 60 100];

v_pdf = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/m));

% Generate an initial population
for i = 1:population_size
    angle = rand*2*pi;
    state(i,:) = [length*rand width*rand random(v_pdf) random(v_pdf)];
    
    % Make sure no particles start in a box
    while(in_box(state(i,1:2), boxes))
        state(i,1:2) = [length*rand width*rand];
    end
end

for i = 1:iterations
    new_state = state;
    
    %Update positions
    new_state(:,1:2) = state(:,1:2) + time_step.*state(:,3:4);
    
    j = state(:,1) > length;
    new_state(j,1) = state(j,1) - length;
    
    j = state(:,1) < 0;
    new_state(j,1) = state(j,1) + length;
    
    j = state(:,2) > width;
    new_state(j,2) = 2*width - state(j,2);
    new_state(j,4) = -state(j,4);
    
    j = state(:,2) < 0;
    new_state(j,2) = -state(j,2);
    new_state(j,4) = -state(j,4);
    
    % Look for particles that have "entered" a box and move them to
    % where they should be. Note that boxes should not overlap for this
    % to work
    for j=1:population_size
        box_num = in_box(new_state(j,1:2), boxes);
        if(box_num == 0) continue;
        end
        
        if(state(j,3) > 0 && state(j,1) > boxes(box_num, 1))
            new_state(j,1) = 2*boxes(box_num, 1) - state(j,1);
            new_state(j,3) = -state(j,3);
        end
        if(state(j,3) < 0 && state(j,1) < boxes(box_num,2))
            new_state(j,2) = 2*boxes(box_num,2) - state(j,1);
            new_state(j,3) = -state(j,3);
        end
        
        if(state(j,4) > 0 && state(j,2) > boxes(box_num, 3))
            new_state(j,1) = 2*boxes(box_num, 3) - state(j,2);
            new_state(j,4) = -state(j,4);
        end
        if(state(j,4) < 0 && state(j,2) < boxes(box_num,4))
            new_state(j,2) = 2*boxes(box_num,4) - state(j,2);
            new_state(j,4) = -state(j,4);
        end
    end
    
    
    % Scatter particles
    j = rand(population_size, 1) < p_scat;
    new_state(j,3:4) = random(v_pdf, [sum(j),2]);
    
    state = new_state;
    
    % Record the temperature
    temperature(i) = (sum(state(:,3).^2) + sum(state(:,4).^2))*m/k/2/population_size;
    
    % Record positions for subset of particles that will be graphed
    % if show_movie = false
    for j=1:plot_population
        trajectories(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    
    % Update the movie every 5 iterations
    if show_movie && mod(i,5) == 0
        figure(1);
        subplot(2,1,1);
        hold off;
        plot(state(1:plot_population,1)./1e-9, state(1:plot_population,2)./1e-9, 'o');
        hold on;
        
        % Plot the boxes
        for j=1:size(boxes,1)
           plot([boxes(j, 1) boxes(j, 1) boxes(j, 2) boxes(j, 2) boxes(j, 1)]./1e-9,...
               [boxes(j, 3) boxes(j, 4) boxes(j, 4) boxes(j, 3) boxes(j, 3)]./1e-9, 'k-');
        end
        
        axis([0 length/1e-9 0 width/1e-9]);
        title(sprintf('Trajectories for %d of %d Electrons with Fixed Velocity',...
        plot_population, population_size));
        xlabel('x (nm)');
        ylabel('y (nm)');
        if i > 1
            subplot(2,1,2);
            hold off;
            plot(time_step*(0:i-1), temperature(1:i));
            axis([0 time_step*iterations min(temperature(1:i))*0.98 max(temperature)*1.02]);
            title('Semiconductor Temperature');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
        end
        pause(0.05);
    end
end

% Show trajectories after the movie is over
figure(1);
subplot(2,1,1);
title(sprintf('Electron Trajectories for %d of %d Electrons with Fixed Velocity',...
    plot_population, population_size));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 length/1e-9 0 width/1e-9]);
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
    subplot(2,1,2);
    hold off;
    plot(time_step*(0:iterations-1), temperature);
    axis([0 time_step*iterations min(temperature)*0.98 max(temperature)*1.02]);
    title('Semiconductor Temperature');
    xlabel('Time (s)');
    ylabel('Temperature (K)');
end

% Store the positions and velocities of the electrons