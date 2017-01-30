%% Part 2: Collisions with Mean Free Path

% Calculate the scattering probability in one time step
p_scat = 1 - exp(-time_step/0.2e-12);

% Each row corresponds to an electron with the position and velocities
% [x y vx vy]
state = zeros(population_size, 4);
trajectories = zeros(iterations, plot_population*2);
temperature = zeros(iterations,1);

% The distribution of velocities in x and y is Gaussian, with a 
% standard deviation of sqrt(k*T/m). This results in an overall
% Maxwell-Boltzmann velocity distribution at temperature T
v_pdf = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/m));

% Generate an initial population
for i = 1:population_size
    angle = rand*2*pi;
    state(i,:) = [length*rand height*rand random(v_pdf) random(v_pdf)];
end

%%
% The average velocity should be calculated to be correct.

avg = sqrt(sum(state(:,3).^2)/population_size + sum(state(:,4).^2)/population_size)

%%
% This returns a velocity of about 18.7 km/s, which is correct.
% This varies a little bit, since the initial velocities are random with a 
% MB distribution.

for i = 1:iterations
    %Update positions
    state(:,1:2) = state(:,1:2) + time_step.*state(:,3:4);
    
    j = state(:,1) > length;
    state(j,1) = state(j,1) - length;
    
    j = state(:,1) < 0;
    state(j,1) = state(j,1) + length;
    
    j = state(:,2) > height;
    state(j,2) = 2*height - state(j,2);
    state(j,4) = -state(j,4);
    
    j = state(:,2) < 0;
    state(j,2) = -state(j,2);
    state(j,4) = -state(j,4);
    
    % Scatter particles
    j = rand(population_size, 1) < p_scat;
    state(j,3:4) = random(v_pdf, [sum(j),2]);
    
    % Record the temperature
    temperature(i) = (sum(state(:,3).^2) + sum(state(:,4).^2))*m/k/2/population_size;
    
    % Record positions for subset of particles that will be graphed
    % if show_movie = 0
    for j=1:plot_population
        trajectories(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    
    % Update the movie every 5 iterations
    if show_movie && mod(i,5) == 0
        figure(2);
        subplot(3,1,1);
        hold off;
        plot(state(1:plot_population,1)./1e-9, state(1:plot_population,2)./1e-9, 'o');
        axis([0 length/1e-9 0 height/1e-9]);
        title(sprintf('Trajectories for %d of %d Electrons (Part 2)',...
        plot_population, population_size));
        xlabel('x (nm)');
        ylabel('y (nm)');
        if i > 1
            subplot(3,1,2);
            hold off;
            plot(time_step*(0:i-1), temperature(1:i));
            axis([0 time_step*iterations min(temperature)*0.98 max(temperature)*1.02]);
            title('Semiconductor Temperature');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
        end
        
        % Show histogram of speeds
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
figure(2);
subplot(3,1,1);
title(sprintf('Trajectories for %d of %d Electrons (Part 2)',...
    plot_population, population_size));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 length/1e-9 0 height/1e-9]);
hold on;
for i=1:plot_population
    plot(trajectories(:,i*2)./1e-9, trajectories(:,i*2+1)./1e-9, '.');
end

% Show temperature plot over time
if(~show_movie)
    subplot(3,1,2);
    hold off;
    plot(time_step*(0:iterations-1), temperature);
    axis([0 time_step*iterations min(temperature)*0.98 max(temperature)*1.02]);
    title('Semiconductor Temperature');
    xlabel('Time (s)');
    ylabel('Temperature (K)');
end

% Show speed histogram
subplot(3,1,3);
v = sqrt(state(:,3).^2 + state(:,4).^2);
title('Histogram of Electron Speeds');
histogram(v);
xlabel('Speed (m/s)');
ylabel('Number of particles');