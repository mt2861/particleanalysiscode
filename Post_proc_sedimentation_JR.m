addpath('./lib')
clear all
close all

set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',22)

lw = 2; % default linewidth
ms = 8; % default markersize
%% Options
image_sequence = 0;
step = 10;
flip_axis=1;
moving_grid = 1;
if image_sequence
    moving_grid = 0;
end
if ~moving_grid
    pixel_size = 5.89/4;
    window = [1216 1936] * pixel_size;
end

save_paraview = 0;

histogram_height = 1;
histogram_height_front = 1;
show_histogram_height = 1;
save_figure_height_histogram =0;
save_movie_histogram = 1;
save_avi = 0;
save_png = 1;

analyze_fourier = 0;
show_movie_fourier = 1;
show_particles = 1;
save_movie_fourier = 1;

analyze_front_density = 1;
shifted_density = 0; % work in progress: bins the particles by both row and 
%and column and then shifts the columns so that all column-maxima occur in
%the same row. Then creates a histogram by averaging across the rows. 
plot_density = 0;
save_movie_density = 0;
front_speed = 0; %not tested

resolution = 0; % 0: singleblob, 12: 12 blobs, etc...
%% declare the remaining needed variables
if image_sequence
    imseq_name = '../20180725/threshold/50';
    frame_rate = 10/step;
    dt = 1/frame_rate;
    pixel_size = 5.89/4;
    % No height data in image sequence
    histogram_height = 0;
    a = 2.1;
else
    N = 30000;
    % particle radius 
    a = 2.1; 
    % Realization number
    Real = '1'; 
    % Gravitational height over a
    hg_over_a = '0_027394'; 
    % TEMPORARY: Fx = 6*pi*a*eta*Vel, In the future we will just use angle 
    % theta on gravity
    Vel = '5a'; 
    string_more ='';
    string_more2 = '';
    initial_x_length = 645;
    initial_y_length = 1611;
    tilt_angle = '1_396';
    
    keyword = ''; %leave blank if no keyword.
    if length(keyword) ==0
    else
        keyword = [keyword '_'];
    end
end
%% Extract parameters from input file
% creates a structure array whos fields are named dynamically according to
% the input file
if image_sequence
    imseq = Imseq(imseq_name, step, flip_axis);
    box_x = imseq.dimensions(1)*pixel_size;
    box_y = imseq.dimensions(2)*pixel_size;
    box_z = 0;
else
    filename_input = [keyword tilt_angle '_angle_' num2str(N) '_N_' num2str(initial_x_length) '_Lx_' num2str(initial_y_length) '_Ly.inputfile'];
    dat = readInputfile( filename_input);

    blob_radius = dat.blob_radius;
    kT = dat.kT;
    eps_w = dat.repulsion_strength_wall;
    eps_p = dat.repulsion_strength;
    str_eps_p = num2str_custom(round(eps_p/kT));
    str_eps_w = num2str_custom(round(eps_w/kT));
    debye_w = dat.debye_length_wall;
    debye_p = dat.debye_length;
    str_debye_p = num2str_custom(debye_p/a);
    str_debye_w = num2str_custom(debye_w/a);
    dt = dat.dt;
    num_steps = dat.n_steps; %number of steps
    step_save_frequency = dat.n_save; %save frequency
    eta = dat.eta;
    mg = dat.g;
    if moving_grid
        box_x = dat.lx;
        box_y = dat.ly;
        box_z = dat.lz;
    else
        box_x = 1216*pixel_size;
        box_y = 1936*pixel_size;
    end

    area_fraction = (N*pi*a^2)/(initial_x_length*box_y);
    str_area_fraction = num2str_custom(area_fraction);
    grav_height = kT/mg/a;
    str_grav_height = num2str_custom(grav_height);
    %tilt_angle = dat.tilt_angle;
    %str_tilt_angle = strrep(num2str(tilt_angle),'.','_');

    % Test blob radius when resolution = 0
    if resolution == 0
       if blob_radius ~= a
           disp('Error blob radius does not correspond to inputfile')
       end
    end

    % Test gravitational height
    hg_input_over_a = kT/mg/a;
    if hg_input_over_a ~= str2num_custom(hg_over_a)
       str2num_custom(hg_over_a)
    end

    num_saved_steps = (num_steps/step_save_frequency);

    dir_to_take_data = '/stoch/JRCruise/RigidMultiblobsWall/multi_bodies/examples/sedimenting_spheres/saves';
    folder = [keyword 'N_' num2str(N) '_area_fraction_' str_area_fraction '_grav_height_' str_grav_height '_Lx_' num2str(initial_x_length) '_Ly_' num2str(box_y) '_tilt_angle_' tilt_angle];
end
%% Defines the physcial time vector and the frames to import
if image_sequence
    saved_steps_vec = 1:imseq.no_frames;
    time_vec = (saved_steps_vec-1)*dt;
    num_time_samples = length(saved_steps_vec);
else
    Still_running = 0;
    if save_paraview == 0
        if Still_running==1
         jump_saves_import = 1;
        else
            jump_saves_import = 1;
        end
    else
        jump_saves_import = 2;
    end
    % normally like arange(Nsaves)
    saved_steps_vec = 1:jump_saves_import:num_saved_steps;
    time_vec = (saved_steps_vec-1)*step_save_frequency*dt; 
    % number of time samples
    num_time_samples = length(saved_steps_vec);
end
%% Open data files
dir_to_save_data = '/stoch/JRCruise/RigidMultiblobsWall/multi_bodies/examples/sedimenting_spheres/saves/';
if image_sequence
else
    file_blobs_tot = [dir_to_take_data '/' folder '/' folder '.config'];
    if exist(file_blobs_tot, 'file')
        data_pos_orient = readConfig(file_blobs_tot, num_saved_steps, N);
        if Still_running == 0
                files_to_move_list = dir( fullfile(dir_to_take_data,[folder '.*']) );
                files_to_move = {files_to_move_list.name};
            for n=1:length(files_to_move)
                movefile([dir_to_take_data files_to_move{n}], dir_to_save_data);
            end
        end
    end
end
%% Check that sizes correspond
if image_sequence
else
    Nsaves_test = size(data_pos_orient,1)/N;
    if Nsaves_test ~= num_saved_steps
        Nsaves_test
        num_saved_steps
        pause
    end
    Ntest = size(data_pos_orient,1)/num_saved_steps 
    if Ntest~=N
      Ntest
      N
      pause
    end
end
%% Define destination for paraview files
my_path = pwd;
if save_paraview == 1
    cd(dir_to_save_data)
    dir_paraview = ['./csv_files']
    mkdir(dir_paraview)
    cd(my_path)
end
%% Declare variables
% Positions of the CoM
if image_sequence
else
    pos_COMs = zeros(num_time_samples,N,3);
end
%% Import data    
% Starts the time loop
% collapse into 3d array
if image_sequence
else
    t=0;
    for n = saved_steps_vec
        t=t+1;
          for l = 1:N
            pos_COMs(t,l,1:3) = data_pos_orient(N*(n-1)+l,1:3);
          end  
    end
end
pos_COMs(:,:,1) = pos_COMs(:,:,1) - min(min(pos_COMs(:,:,1)));
%% Analyze density profile along translational direction
if front_speed
    peak_positions = [];
end
if analyze_front_density ==1
    if image_sequence
        saves_to_analyze = saved_steps_vec;
    else
        list_part = 1:N;
        hmin = min(min(pos_COMs(:,:,3)));
        hmax = max(max(pos_COMs(:,:,3)));
        jump_front_density = 1;
        save_1 = 1;
        save_end = num_time_samples;   
        saves_to_analyze = save_1:jump_front_density:save_end;
    end  
    % Number of frames/saved states
    Nanalyzed = length(saves_to_analyze);
    % Bin
    % if moving grid dynamically choose grid and gridsize
    if moving_grid
        max_x = max(max(pos_COMs(save_1:save_end,:,1)));  
        min_x = min(min(pos_COMs(save_1:save_end,:,1)));
        dx_density = 5*a;
        max_y = max(max(pos_COMs(save_1:save_end,:,2)));
        min_y = min(min(pos_COMs(save_1:save_end,:,2)));
        dy_density = 15*a;
    % else choose a grid aligned with pixel grid
    else
        max_x = window(1);
        min_x = 0;
        dx_density = pixel_size;
    end
    % Set required variables from image sequence properties
    if image_sequence
        max_x = imseq.dimensions(1)*pixel_size;
        min_x = 0;
        dx_density = pixel_size;
    end
    % Common variables
    if image_sequence
        x_density = (0:imseq.dimensions(1)-1)*pixel_size;
        density_x_time = zeros(Nanalyzed,imseq.dimensions(1));
    else
        Nx_density = floor((max_x-min_x)/dx_density);
        x_density = linspace(min_x,max_x,Nx_density);
        dx_density = x_density(2)-x_density(1);
        density_x_time = zeros(Nanalyzed,Nx_density);
        if shifted_density
            Ny_density = floor((max_y-min_y)/dy_density);
            y_density = linspace(min_y,max_y,Ny_density);
            dy_density = y_density(2)-y_density(1);
            shifted_number_grid = struct;
        end
    end
    % set variable density_x_time
%     wait = waitbar(0)
    if image_sequence
        for i=saves_to_analyze
            imseq_i = imseq.get(i);
            thresholded = imseq_i == 0;
            density_x_time(i, 1:end) = sum(imseq_i, 2);
            waitbar(i/imseq.no_frames, wait,sprintf('frame %d of %d', i, imseq.no_frames));
        end
        disp(density_x_time(1,1:end))
    else
        t = 0;
        for i=saves_to_analyze
            t=t+1;
            indices_x = floor(pos_COMs(i,:,1)/dx_density)+1;
            for nx=1:Nx_density
                density_x_time(t,nx) = numel(list_part(indices_x == nx));
            end
            if front_speed
                [dumb,peak_position] = max(density_x_time(t,:))
                peak_positions = [peak_positions peak_position*dx_density+dx_density/2];
            end
            if shifted_density
                number_grid = zeros(Nx_density,Ny_density);
                for j = list_part
                    index = pos_COMs(i,j,1:2);
                    grid_y = floor(index(2)/dy_density)+1+abs(floor(min_y/dy_density));
                    grid_x = floor(index(1)/dx_density)+1;
                    number_grid(grid_x,grid_y) = number_grid(grid_x,grid_y) + 1;
                end
                pos_max = zeros(1,Ny_density);
                for l = 1:Ny_density
                    cont = true;
                    m = grid_x;
                    while cont
                        if (m+1 > size(number_grid,1) || number_grid(m+1,l) == 0)
                            cont = false;
                            pos_max(l) = m;
                        end
                        m = m+1;
                    end
                end
                II = min(pos_max);
                shift = max(pos_max)-min(pos_max);
                shftd_expd_grid = [number_grid ; zeros(shift,Ny_density)];
    %             shftd_expd_grid = [zeros(shift,Ny_density) ; number_grid];
                for k = 1:Ny_density
                    shftd_expd_grid(:,k) = circshift(shftd_expd_grid(:,k),pos_max(k)-II);
                end
                shifted_histograms.(['save' num2str(t)]) = struct;
                shifted_histograms.(['save' num2str(t)]).grid = shftd_expd_grid;
                shifted_histograms.(['save' num2str(t)]).averaged_data = sum(shftd_expd_grid,2);
                shifted_histograms.(['save' num2str(t)]).shift = shift;
            end
        end
    end

    % Plot
    if plot_density == 1
        hfig = figure;
        set(hfig,'position',[500 500  1500 800])
        % Run through frames
        t = 0;    
        for i=saves_to_analyze
            t=t+1;    
            % Settings
            box on; subplot(1,3,1); box on
            % Plot COMs
            if image_sequence
                axis_params = [0 box_y/a min_x/a max_x/a];
                frame = imseq.get(i);
                plotImage(i, t, frame, pixel_size, a, box_x, box_y, ms, axis_params, time_vec)
            else
                if ~moving_grid
                    axis_params = [0 window(2)/a min_x/a max_x/a];
                else
                    axis_params = [0 box_y/a min_x/a max_x/a];
                end
                plotCOM(i, t, pos_COMs, a, box_x, box_y, ms, axis_params, ...
                    hmin, hmax, time_vec)
            end
            axis equal
            axis(axis_params)
            xlabel('$y/a$','fontsize',30)
            ylabel('$x/a$','fontsize',30) 
            set(gca,'yminortick','on')
            set(gca,'xminortick','on')
            set(gca,'ticklength',3*get(gca,'ticklength'))
            box on
            title(['$t = $' num2str(time_vec(i),'%5.2f') '$ s$'],'fontsize',30)
            % Settings
            subplot(1,3,2:3); box on
            % Plot density
            plot((x_density + dx_density/2)/a,density_x_time(t,:),'-k','linewidth',1.5)    
            xlim([min_x/a max_x/a])
            ylim([0 max(max(density_x_time))])     
            % Agnostic settings
            ylabel('$\rho(x,t)$','fontsize',30)     
            xlabel('$x/a$','fontsize',30)
            set(gca,'yminortick','on')
            set(gca,'xminortick','on')
            set(gca,'ticklength',3*get(gca,'ticklength'))
            set(hfig,'color','w')
            set(gca,'layer','top')
            pause(0.0001)
            % Save frame
            if save_movie_density == 1
                mov(t) = getframe(hfig);
                if t==1
                    size_mov1 = size(mov(t).cdata,1);
                    size_mov2 = size(mov(t).cdata,2);
                end
                mov(t).cdata =  mov(t).cdata(1:size_mov1-10,1:size_mov2-10,:);
            end
            % Clear frame if not the last one
            if t<Nanalyzed; clf; end
        end
    end
    % Save movie
    if save_movie_density == 1
        filename_movie = 'Density_along_x_';
        saveMovie(dir_to_save_data, filename_movie, mov )
    end
    % Save variables
    to_save = [[-1 x_density]; [time_vec(saves_to_analyze)', density_x_time]];
    filename = ['Density_x_along_time_dx_' num2str(dx_density/a) 'a'];
    save([dir_to_save_data '/' filename '.txt'],'to_save','-ascii')  
end
if front_speed 
    delta_t = time_vec(2)-time_vec(1)
    peak_velocities(1) = 0;
    for i = 2:length(time_vec)
        peak_velocities(i) = (peak_positions(i)-peak_positions(i-1))/delta_t
    end
    filename = [dir_to_save_data '/' folder '/front_pos.txt'];
    fid = fopen(filename,'w');
    fprintf(fid,'%s %s %s\n',['time_steps' 'position' 'velocity']);
    fprintf(fid,'%f %f %f\n',[time_vec peak_positions peak_velocities]);
end
%% Histogram height + movie from top
if histogram_height == 1
    t=0; 
    % Max/min position over all time
    hmin = min(min(pos_COMs(:,:,3)));
    hmax = max(max(pos_COMs(:,:,3)));
    % set length of box depending on periodicity
    if box_x == 0
        xmin =  min(min(pos_COMs(:,:,1)));
        xmax =  max(max(pos_COMs(:,:,1)));
    else
        xmin =  0;
        xmax =  box_x;
    end
    if box_y == 0
        ymin =  min(min(pos_COMs(:,:,2)));
        ymax =  max(max(pos_COMs(:,:,2)));
    else
        ymin =  0;
        ymax =  box_y;
    end 
    % set bin width
    bw = 0.1*a;
    Nbin = ceil((hmax-hmin)/bw)
    mean_pdf_height = zeros(1,Nbin);
    mean_pdf_height_front = zeros(1,Nbin);
    pdf_height_time = zeros(num_time_samples,Nbin);
    pdf_height_time_front = zeros(num_time_samples,Nbin);
    mean_height = zeros(num_time_samples,1);
    mean_height_front = zeros(num_time_samples,1);
    
    if show_histogram_height ==1
     hfig=figure;
     set(hfig,'position',[500 500  1200 1000])
    end
    
    for n = saved_steps_vec
        t=t+1 
        
        % NOT RELEVANT FOR MOVIE INPUT        
        % Histogram of particle in the front
        if histogram_height_front ==1
            % Defines the front according to the percentage of particles in the
            % region
            percentage_front = 0.2;
            num_ps_in_front = floor(percentage_front*N);
            tol_front =  1e-3; % has to be within 1e-3 of 0.5
            step_move_front = 0.1;
            [x_curr, indicies_ps_in_front] = ...
                getFront(pos_COMs(t,:,1), num_ps_in_front, tol_front, step_move_front );
            % histogram of particles in front
            [hist_height_front,edges_front] = histcounts( ...
                pos_COMs(t,indicies_ps_in_front,3), ... 
                'BinWidth', bw, ...
                'BinLimits', [hmin,hmax] ...
                );
            [pdf_height_front, xhist_front] = hist2pdf(hist_height_front, edges_front);
            % store and update
            pdf_height_time_front(t,:) = pdf_height_front;
            mean_pdf_height_front = mean_pdf_height_front + pdf_height_front;
            mean_height_front(t) = trapz(xhist_front,xhist_front.*pdf_height_front );     
        end

       % Height Histrogram for all particles
        [hist_height,edges] = histcounts( ...
            pos_COMs(t,:,3), ...
            'BinWidth',bw, ...
            'BinLimits',[hmin,hmax] ...
            );
        [pdf_height, xhist] = hist2pdf(hist_height, edges);
        
        thresholded_heights = hist_height;
        thresholded_heights(thresholded_heights<7) = 0;
        thresholded_heights = thresholded_heights(1:find(thresholded_heights,1,'last'));
        thresholded_max_height = length(thresholded_heights)*bw+bw/2+hmin;
        
        
        % save initial distribution
        if t==1
            xhist_ini  = xhist;
            pdf_height_ini = pdf_height;
        end
          
        pdf_height_time(t,:) = pdf_height;
        mean_pdf_height = mean_pdf_height + pdf_height;        
        mean_height(t) = trapz(xhist,xhist.*pdf_height );
        
        if show_histogram_height ==1
            % Show particles
            subplot(4,1,1:2); box on; hold on;
            mean_x = mean(mod(pos_COMs(t,:,1),box_x));
            max_x = max(pos_COMs(t,:,1));
            axis_params = [ymin/a ymax/a max_x/a-280 max_x/a+20 ];
            plotCOM( n, t, pos_COMs, a, box_x, box_y, ms, axis_params, hmin, thresholded_max_height, time_vec); 
            axis equal
            axis(axis_params)
%             xlabel('$y/a$','fontsize',30)
            ylabel('$x/a$','fontsize',30) 
            set(gca,'yminortick','on')
            set(gca,'xminortick','on')
            set(gca,'ticklength',3*get(gca,'ticklength'))
            box on
            title(['$t = $' num2str(time_vec(t),'%5.2f') '$ s$'],'fontsize',30)
            if  histogram_height_front ==1
                % Show front dividing line 
                plot(linspace(ymin,ymax,100)/a, x_curr*ones(1,100)/a, '--r', 'linewidth', 3);
            end  
            % Show histogram
            subplot(4,1,3); box on; hold on
            plot(xhist/a, pdf_height*a, '-k', 'linewidth',2)    
            plot(xhist_ini/a, pdf_height_ini*a, '-', 'color', [0.7 0.7 0.7], 'linewidth',2)      
            plot(mean_height(t)*ones(10)/a, linspace(0,1,10), '--k', 'linewidth',2)
            
            % Show front height pdf
            if  histogram_height_front ==1
                plot(xhist_front/a, pdf_height_front*a, '-r', 'linewidth',2)
                plot(mean_height_front(t)*ones(10)/a, linspace(0,1,10), '--r', 'linewidth',2)
            end
            
            xlabel('$h/a$','fontsize',24)
            ylabel('$P(h)\times a$','fontsize',24)
%             axis([1 hmax/a 1e-3 3])
%             axis([0 2 0 10])
            axis([0 thresholded_max_height/a 0 20])
            set(hfig,'color','w')
            set(gca,'yscale','log')
            set(gca,'layer','top')
            set(gca,'yminortick','on')
            set(gca,'xminortick','on')
            set(gca,'ticklength',3*get(gca,'ticklength'))
            
            %show packing fraction
            subplot(4,1,4);
            box on 
            hold on
%             plot((x_density + dx_density/2)/a,density_x_time(t,:),'-k','linewidth',1.5) %number density
            plot((x_density + dx_density/2)/a,density_x_time(t,:)*pi*a^2/(dx_density*(max_y-min_y)),'-k','linewidth',1.5) %packing fraction
            if shifted_density
                plot(([-shifted_histograms.(['save' num2str(t)]).shift*dx_density:dx_density:-dx_density x_density]+dx_density/2+abs(min_y))/a, shifted_histograms.(['save' num2str(t)]).averaged_data*pi*a^2/(dx_density*(max_y-min_y)),'r','linewidth',1.5)
            end 
            xlim([min_x/a max_x/a])
            ylim([0 max(max(density_x_time))*pi*a^2/(dx_density*(max_y-min_y))])
            % Agnostic settings
            ylabel('$\phi(x,t)$','fontsize',30)     
            xlabel('$x/a$','fontsize',30)
            set(gca,'yminortick','on')
            set(gca,'xminortick','on')
            set(gca,'ticklength',3*get(gca,'ticklength'))
            set(hfig,'color','w')
            set(gca,'layer','top')
            pause(.0001)
        
            if save_movie_histogram == 1
                if save_avi
                    mov_histo(t) = getframe(hfig);
                    if t==1
                        size_mov1 = size(mov_histo(t).cdata,1)
                        size_mov2 = size(mov_histo(t).cdata,2)
                    end
                    mov_histo(t).cdata =  mov_histo(t).cdata(1:size_mov1-10,1:size_mov2-10,:);
                end
                if save_png
                    frame_enumeration_length = length(num2str(max(saved_steps_vec)));
                    frame_number = num2str(t - 1);
                    while length(frame_number) < frame_enumeration_length
                        frame_number = ['0' frame_number];
                    end
                    saveas(hfig,[dir_to_save_data '/' folder '/Histogram_frame_' frame_number '.png']);
                end
            end
            pause(0.0001)
            if t<num_time_samples
                clf
            end
        end
    end
    
    if save_movie_histogram == 1
        if save_avi
            pwd
            filename_movie = 'Histogram_height';
            saveMovie( dir_to_save_data, filename_movie, mov_histo )
        end
        if save_png == 1
            cd([dir_to_save_data folder])
            if exist('film')
                unix('rm -r -f film');
            end 
            mkdir film
            unix(['convert -resize 50% *.png ' folder '.gif']);
            unix('find -name "Histogram_frame*" -exec mv -t film {} +')
            cd('../..');
        end
    end
    
    mean_pdf_height = mean_pdf_height/num_time_samples;
    if  histogram_height_front ==1
        mean_pdf_height_front = mean_pdf_height_front/num_time_samples;
    end
    
    if save_figure_height_histogram == 1
        hfig = figure;
        box on
        subplot(2,1,1)
        hold on
        plot(xhist/a,mean_pdf_height*a,'-k','linewidth',2)
        if  histogram_height_front ==1
            plot(xhist_front/a,mean_pdf_height_front*a,'-r','linewidth',2)
        end
        xlabel('$h/a$','fontsize',24)
        ylabel('$P(h)\times a$','fontsize',24)
        xlim([0 50])
        set(hfig,'color','w')
        set(hfig,'position',[500 500  1200 1000])
        set(gca,'layer','top')
        set(gca,'yminortick','on')
        set(gca,'xminortick','on')
        set(gca,'ticklength',3*get(gca,'ticklength'))


        subplot(2,1,2)
        hold on
        plot(xhist/a,mean_pdf_height*a,'-k','linewidth',2)
        if  histogram_height_front ==1
           plot(xhist_front/a,mean_pdf_height_front*a,'-r','linewidth',2)
        end
        set(gca,'xscale','log','yscale','log')
        xlabel('$h/a$','fontsize',24)
        ylabel('$P(h)\times a$','fontsize',24)
        set(hfig,'color','w')
        set(hfig,'position',[500 500  1200 1000])
        set(gca,'layer','top')
        set(gca,'yminortick','on')
        set(gca,'xminortick','on')
        set(gca,'ticklength',3*get(gca,'ticklength'))
        title(['Averaged over $t = $' num2str(time_vec(saved_steps_vec(end)),'%5.2f') '$ s$'],'fontsize',30)

        filename =[ 'Mean_histogram_height_time'];
        saveas(hfig,[dir_to_save_data '/' filename '.fig'],'fig')
        saveas(hfig,[dir_to_save_data '/' filename '.eps'],'psc2')
    end
     
    to_save = [[-1 xhist]; [time_vec', pdf_height_time]];
    filename = ['PDf_height_along_time'];
    save([dir_to_save_data '/' filename '.txt'],'to_save','-ascii')
    
    if  histogram_height_front ==1
        to_save = [[-1 xhist]; [time_vec', pdf_height_time_front]];
         filename = ['PDf_height_along_time_front'];
        save([dir_to_save_data '/' filename '.txt'],'to_save','-ascii')
    end
end

%% Fourier analysis of the front along the y-direction
if analyze_fourier == 1 
    if image_sequence
    else
        list_part = 1:N;
        hmin = min(min(pos_COMs(:,:,3)));
        hmax = max(max(pos_COMs(:,:,3)));
    end
    % Choose number of modes such that the minimum wavelength is lambda_min
    if image_sequence 
        Ngrid_fft = imseq.dimensions(2);
        lambda_min = box_y * 2/Ngrid_fft;
        num_a = lambda_min/a;
    elseif ~moving_grid
        Ngrid_fft = 1936;
    else
        lambda_min = 5*a;
        Ngrid_fft = floor(box_y/lambda_min)*2;
        if mod(Ngrid_fft,2) == 1
           Ngrid_fft = Ngrid_fft+1; 
        end
    end
    % Define grid in real space
    y_grid = linspace(0,box_y,Ngrid_fft);
    dy = y_grid(2)-y_grid(1);
    % Sampling frequency
    Fs = 1/dy;                    
    dFy = Fs/Ngrid_fft;  
    Fy = -Fs/2:dFy:Fs/2-dFy;
    % Choose how many saves to analyze
    jump_fourier = 1;
%     saves_to_analyze = 1:jump_fourier:floor(num_time_samples)/4;
    saves_to_analyze = 1:jump_fourier:floor(num_time_samples);
    Nanalyzed = length(saves_to_analyze);
    % Declare arrays
    PSD_time = zeros(Nanalyzed,Ngrid_fft/2);
    autocorrelation_time = zeros(Nanalyzed,Ngrid_fft);
    autocorrelation_distance = zeros(Nanalyzed,1);
    Npart_current = zeros(Nanalyzed,1);
    lower_x = zeros(Nanalyzed,1);
    first_moment_front = zeros(Nanalyzed,1);
    heaviest_freq_front = zeros(Nanalyzed,1);
    % Create figure to display movie
    if show_movie_fourier==1
        hfig800=figure(800);
        set(hfig800,'color','w')
        set(hfig800,'position',[500 500 1000 1200]) 
    end    
    % Define parameters to chose front
    if moving_grid; dx = 0.2*dy; else; dx = dy; end
    if image_sequence
        N = sum(sum(imseq.get(1)));
    end
    percentage_particles_front =floor(0.5*N);
    if moving_grid
        tol_particles = floor(0.005*dx/0.7*N);
    else
        tol_particles = floor(0.005*dx/0.7*N);
    end

    % Main Loop    
    t = 0;
    for i=saves_to_analyze 
        t=t+1
        % Set up parameters
        if image_sequence
            % Set x_grid
            N = sum(sum(imseq.get(i)));
            x_grid = (1:imseq.dimensions(1))*pixel_size;
            % Set lower_x
            frame = imseq.get(i);
            dim1 = imseq.dimensions(1);
            if t==1; prev = dim1; else; prev = lower_x(t-1); end
            NpartFunc = @(lower_x) sum(sum(frame(lower_x:dim1,1:end)));
        else
            indices_y = ceil(pos_COMs(i,:,2)/dy);
            if ~moving_grid
                min_x_grid = 0;
                max_x_grid = 1216*pixel_size;
            else
                % Defines grid along x to decide the portion of particles 
                % on which to average along x
                min_x_grid = min(pos_COMs(i,:,1))-dx;
                max_x_grid = max(pos_COMs(i,:,1));
            end
            x_grid = min_x_grid:dx:max_x_grid;
            Ngrid_x = length(x_grid);
            indices_x = ceil((pos_COMs(i,:,1)-min_x_grid)/dx);
            % Find lower_x
            if t==1; prev = 0; else; prev = lower_x(t-1); end
            NpartFunc = @(lower_x) sum(indices_x >= lower_x);
        end
        % Find lower
        [lower_x_, Npart_current_] = findLower(...
            prev, NpartFunc, percentage_particles_front, tol_particles);
        lower_x(t) = lower_x_;
        Npart_current(t) = Npart_current_;
        % Set density part
        if image_sequence
            density_part = frame(lower_x(t):dim1,1:end);
        else
            density_part = zeros(Ngrid_x-lower_x(t)+1, Ngrid_fft);
            for n = 1:Ngrid_fft
                if show_particles
                    density_part(1, n) = sum(indices_y==n);
                else
                    for nx = 1:Ngrid_x-lower_x(t)+1
                        density_part(nx,n) = numel(list_part( indices_y==n & indices_x==(nx+lower_x(t)-1) ));
                    end
                end 
            end
        end
        density_to_analyze = sum(density_part,1);       
        signal_front = (density_to_analyze-mean(density_to_analyze)); 
        % autocorrelation function
        [acf,lags] = xcorr(signal_front,'coeff');
        acf = acf(lags>=0); % Only keep autocorelation for positive lags
        acffilter = sgolayfilt(acf,2,9); % smooth the autocorr
        autocorrelation_time(t,:) = acffilter;        
        ind_start = find(y_grid/a<=30,1,'last');
        ind_end = find(y_grid/a<=100,1,'last');    

        [val,I] = max(acffilter(ind_start:ind_end));
        autocorrelation_distance(t) = y_grid(ind_start+I-1)/a;
        % FFT analysis
        Y_front = fftshift(fft(signal_front)/Ngrid_fft);        
        Fy_half = Fy(Ngrid_fft/2+1:Ngrid_fft);
        Y_front_half = Y_front(Ngrid_fft/2+1:Ngrid_fft);
        PSD_time(t,:) = abs(Y_front_half).^2;
        % sort frequencies by magnitude store in Fsort_front
        [Spectrum_front_sorted,ind_freq]=sort(abs(Y_front_half),'descend');
        Fsort_front = abs(Fy_half(ind_freq));
        
        Normalized_distrib_front = abs(Y_front_half).^2/trapz(abs(Fy_half),abs(Y_front_half).^2);
        first_moment_front(t) = trapz(abs(Fy_half),Normalized_distrib_front.*abs(Fy_half));
        heaviest_freq_front(t) = Fsort_front(1);
        % find the max freq more than .005
        [none, start_freq_ind] = max(Fy_half > a/200);
        disp(a/(Fy_half(start_freq_ind)));
        [Spectrum_front_sorted,ind_freq]=sort(abs(Y_front_half(start_freq_ind:end)),'descend');
        Fsort_front = abs(Fy_half(min(length(Y_front_half),ind_freq+start_freq_ind)));
        %% Plot image
        subplot(4,2,1:2); box on; hold on
        title(['$t = $' num2str(time_vec(i),'%5.2f') '$ s$'],'fontsize',30)
%         axis_params = [y_grid(1)/a y_grid(end)/a 0 x_grid(end)/a];
        axis_params = [y_grid(1)/a y_grid(end)/a x_grid(lower_x(t))/a x_grid(end)/a];
        if show_particles == 1
            if image_sequence
                plotImage(i, t, frame, pixel_size, a, box_x, box_y, ms, axis_params, time_vec)
            else
                plotCOM(i, t, pos_COMs, a, box_x, box_y, ms, axis_params, hmin, hmax, time_vec)
            end
            axis equal; shading interp
            xlim([y_grid(1)/a y_grid(end)/a])
            ylim([x_grid(lower_x(t))/a x_grid(end)/a])
%             ylim([0 x_grid(end)/a])
            ylabel('$x/a$','fontsize',30)
            set(gca,{'xminortick', 'ticklength', 'layer'} ,{'on', 1.5*get(gca,'ticklength'), 'top'})
        else
            pcolor(y_grid/a,x_grid(lower_x(t):end)/a,density_part)
            colormap(flipud(autumn))
            caxis([0 3])                   
        end
        %% Plot front
        subplot(4,2,3:4); box on; hold on  
        plot(y_grid/a,signal_front,'-k','linewidth',lw)
        xlim([y_grid(1)/a y_grid(end)/a])
        lim_signal = floor(tol_particles/2);
        ylim([-lim_signal lim_signal])
        ylabel('$\delta \rho$ front','fontsize',30)
        set(gca,{'yminortick', 'xminortick', 'ticklength', 'layer'} ,{'on', 'on', 3*get(gca,'ticklength'), 'top'})
        %% Plot autocorr
        subplot(4,2,5:6); box on; hold on
        plot(y_grid/a, acffilter,'linewidth',lw) ;      
        plot(autocorrelation_distance(t),val,'sk','linewidth',lw)
        text(autocorrelation_distance(t)*1.2,...
             val*1.5,...
             num2str(autocorrelation_distance(t)),...
             'fontsize',20)
        xlim([y_grid(1)/a y_grid(end)/a])
        ylim([-0.5 1])
        xlabel('$y/a$','fontsize',24)
        ylabel('Autocorr','fontsize',30)
        set(gca,{'yminortick', 'xminortick', 'ticklength', 'layer'} ,{'on', 'on', 3*get(gca,'ticklength'), 'top'})
        %% Plot FFT
        ymax = 1;
        subplot(4,2,7:8); box on; hold on
        plot(Fy_half*a,abs(Y_front_half).^2,'linewidth',lw) ;
        plot(Fsort_front(1)*a, min([Spectrum_front_sorted(1).^2 ymax]),'sk','linewidth',lw)
        text(Fsort_front(1)*a*1.15,...
             min([(Spectrum_front_sorted(1).^2)*0.8 ymax]),...
             num2str(1/(Fsort_front(1)*a)),... 
             'fontsize',20)
        plot(first_moment_front(t)*a*ones(1,100), ...
            linspace(0,max(abs(Y_front_half).^2),100), ...
            '-.k', ...
            'linewidth',lw)         
        xlim([min(Fy_half*a) max(Fy_half*a)])
%         ylim([0 0.3])
        ylim([0 ymax])
        xlabel('$a/\lambda$','fontsize',24)
        ylabel('$|\delta \hat{\rho}|^2$','fontsize',30)
        set(gca,{'yminortick', 'xminortick', 'ticklength', 'layer'} ,{'on', 'on', 3*get(gca,'ticklength'), 'top'})
        %% Store frame
        pause(0.00001)
        if save_movie_fourier == 1
            mov_fourier(t) = getframe(hfig800);
            if t==1
                 size_mov_fourier_1 = size(mov_fourier(t).cdata,1);
                 size_mov_fourier_2 = size(mov_fourier(t).cdata,2);
            end
            mov_fourier(t).cdata =  mov_fourier(t).cdata(1:size_mov_fourier_1-10,1:size_mov_fourier_2-10,:);
        end
        if t<Nanalyzed
           clf 
        end
    end
    
    if save_movie_fourier == 1
      filename_movie =[ 'Fourier_analysis_Ngrid_' num2str(Ngrid_fft) ...
                        '_lambda_min_' num2str(lambda_min/a) 'a' ...
                        '_percentage_front_' num2str(percentage_particles_front/N) ...
                        '_dy_over_dx_' num2str(dy/dx)];
       writerObj2 = VideoWriter([dir_to_save_data '/' filename_movie '.avi']);
       writerObj2.Quality = 75;
       writerObj2.FrameRate = 5;
       open(writerObj2)
        for t = 1:length(mov_fourier)
               if size(mov_fourier(t).cdata,1)>0 &&  size(mov_fourier(t).cdata,2)>0
                   writeVideo(writerObj2,mov_fourier(t))
               end
        end
        close(writerObj2)
    end
    
    filename = ['PSD_time_Ngrid_' num2str(Ngrid_fft) ...
                '_lambda_min_' num2str(lambda_min/a) 'a' ...
                '_percentage_front_' num2str(percentage_particles_front/N) ...
                '_dy_over_dx_' num2str(dy/dx)];
            
    to_save = [0 1./(Fy_half*a); time_vec(saves_to_analyze)', PSD_time];
    save([dir_to_save_data '/' filename '.txt'],'to_save','-ascii')
    
    filename = ['Autocorrelation_time_Ngrid_' num2str(Ngrid_fft) ...
                '_lambda_min_' num2str(lambda_min/a) 'a' ...
                '_percentage_front_' num2str(percentage_particles_front/N) ...
                '_dy_over_dx_' num2str(dy/dx)];
            
    to_save = [0 y_grid/a; time_vec(saves_to_analyze)', autocorrelation_time];
    save([dir_to_save_data '/' filename '.txt'],'to_save','-ascii')
    
    filename = ['Autocorrelation_distance_Ngrid_' num2str(Ngrid_fft) ...
                '_lambda_min_' num2str(lambda_min/a) 'a' ...
                '_percentage_front_' num2str(percentage_particles_front/N) ...
                '_dy_over_dx_' num2str(dy/dx)];
            
    to_save = [time_vec(saves_to_analyze)', autocorrelation_distance];
    save([dir_to_save_data '/' filename '.txt'],'to_save','-ascii')
end
% Post_proc_sedimentation_JR.m
% Displaying Post_proc_sedimentation_JR.m.
