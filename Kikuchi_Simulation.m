% Version 1.2.5 - 27 Oct 2015 - Last edits by Adam D. Herron
% Most recently altered band thickness for spherical patterns
classdef Kikuchi_Simulation < handle
    % Kikuchi_NORM is a class for handling diffraction or structure factor
    % data and generating Kikuchi Patterns from them. Data can be read from a data
    % file or imported manually.
    properties
        reciprocal_points
        reciprocal_spacing
        N_atoms
        intensity
        lambda
        filedata
        screen_X
        screen_Y
        screen_Z
        screen_V
        complex
        RE_IM
        zone
        SAED_Intensity
    end
    
    methods
        
        function this = Kikuchi_Simulation(lambda)
            if nargin < 1 || isempty(lambda)
                this.lambda = 0.0251;
            else
                this.lambda = lambda;
            end
        end
        
        function Readfile(this,fname,reciprocal_spacing,N_atoms) % used to import data from an "saed.txt" file
            
            if nargin < 2 || isempty(fname) || ~exist(fname,'file')
                [filename, pathname] = uigetfile( ...
                   {'*.txt', 'Text file (*.txt)';
                    '*.*',  'All Files (*.*)'}, ...
                    'Pick a file');
                if filename==0
                    return
                end
                fname=[pathname,filename];
            end
            
            this.reciprocal_points = [];
            this.intensity = [];
            this.reciprocal_spacing = [];
            this.filedata=importdata(fname,' ',4);
            this.reciprocal_points = this.filedata.data(:,2:4); % read in reciprocal sample lattice
            if nargin < 3 || isempty(reciprocal_spacing)
                this.reciprocal_spacing = ... % assume that the distance between the first two reciprocal_points is representitive of the reciprocal spacing of the sample lattice
                    norm(this.filedata.data(1,2:4)-this.filedata.data(2,2:4)); % This is a bad assumption to make, but the calculation here ins't critical and is only used as a default value if the user doesn't specify one
            else
                this.reciprocal_spacing = reciprocal_spacing;
            end
            if nargin < 4 || isempty(N_atoms)
                this.N_atoms = 108000; % set 108000 as the default number of atoms in a simulation
            else
                this.N_atoms = N_atoms;
            end
            if length(this.filedata.data(1,:)) == 6
                this.complex = false;
                this.intensity(:,1) = real(this.filedata.data(:,6)); % read in intensity values
                % this.intensity(:,2) = log10(this.intensity(:,1)); % calculate log base 10 of intensities and sotre in 2nd column
            elseif length(this.filedata.data(1,:)) == 7
                this.complex = true;
                this.RE_IM(:,1) = this.filedata.data(:,6); % read in values for the real part of structure factor
                this.RE_IM(:,2) = this.filedata.data(:,7); % read in values for the imaginary part of structure factor
            else
                error('"Data file is not recognizable by Kikuchi_Simulation class (too many or too few columns in the data series)"');
            end
        end
        
        function Inherit_Data(this,point_array,intensity_array,reciprocal_spacing,N_atoms) % option for importing data without using .Readfile
            this.reciprocal_points = [];
            this.intensity = [];
            this.reciprocal_spacing = [];
            if length(point_array(1,:))~=3
                error('point_array must be (n x 3)');
            end
            if length(point_array(:,1))~=length(intensity_array(:,1))
                error('point_array and intensity_array must have the same number of rows');
            end
            this.reciprocal_points = point_array;
            if nargin < 4 || isempty(reciprocal_spacing)
                this.reciprocal_spacing = ... % assume that the distance between the first two reciprocal_points is representitive of the reciprocal spacing of the sample lattice
                norm(this.reciprocal_points(1,1:3)-this.reciprocal_points(2,1:3));
            else
                this.reciprocal_spacing = reciprocal_spacing;
            end
            if nargin < 5 || isempty(N_atoms)
                this.N_atoms = 32000;
            else
                this.N_atoms = N_atoms;
            end
            if length(intensity_array(1,:)) == 1
                this.complex = false;
                this.intensity(:,1) = real(intensity_array);
%                 this.intensity(:,2) = log10(intensity_array);
            elseif length(intensity_array(1,:)) == 2
                this.complex = true;
                this.RE_IM = intensity_array;
            else
                error('intensity_array must be an n x 1 array (for intensity values) or an n x 2 array (for structure factor real and imaginary parts)');
            end
        end

        function Kikuchi_Calculation(this,zone,dpi,threshold_magnitude,screendist,complex)
            if nargin < 2 || isempty(zone)
                zone = [0,0,0];
            end
            if length(zone(1,:)) ~= 3
                error('error: zone must be 1x3');
            end
            if nargin < 3 || isempty(dpi)
                dpi = 256;
            end
            if nargin < 4 || isempty(threshold_magnitude)
                threshold_magnitude = 3;
            end
            if nargin < 5 || isempty(screendist)
                screendist = 3; %in
            end
            if nargin < 6 || isempty(complex)
                complex = this.complex;
            end
            resolution_width = round(dpi * 6);
            [X,Y,Z] = this.Screen_Vectors(resolution_width,zone,screendist);
            ilambda = 1/this.lambda;
            K = this.Threshold_Data(threshold_magnitude,complex);
            K_mag = sqrt(K(:,1).^2 + K(:,2).^2 + K(:,3).^2);
            
            % dK = this.reciprocal_spacing; % Actual dK
            dK = 4*ilambda/sqrt(((dpi^2)*(screendist^2)+1)); % dK_eff (lower limit on dK for specified dpi)
            % dK = 4*ilambda*((1+((R_mag.^2)*(this.lambda^2)/4)).^(-3/2))/(screendist*dpi); %Alternate dK_eff
            
            if dK < this.reciprocal_spacing
                dK = this.reciprocal_spacing; % minimum bandwidth for a single sample point is 1-2 pixels
            end
            
            if zone(1) == 0 && zone(2) == 0 && zone(3) == 0
                dK = 2*dK;
            end
            
            this.zone = zone;
            % upper_bound = sqrt(ilambda^2 + (R_mag*dK/2)); % Actual dK
            % lower_bound = sqrt(ilambda^2 - (R_mag*dK/2)); % Actual dK
            upper_bound = sqrt((ilambda^2)+K_mag.*dK/2); % Req dK
            % upper_bound = sqrt((ilambda^2)-K_mag); % All the way to theta == 0
            lower_bound = sqrt((ilambda^2)-K_mag.*dK/2); % Req dK
            if complex
                SF1 = zeros(size(X));
                SF2 = zeros(size(X));
            else
                V = zeros(size(X));
            end
            tic
            W = waitbar(0,'Please Wait...','Name','Classic Kikuchi Calculation');
            if complex
                for K_index = 1:length(K(:,1))
                    dist = sqrt((X-K(K_index,1)).^2+(Y-K(K_index,2)).^2+(Z-K(K_index,3)).^2);
                    SF1_contribution = K(K_index,4)*ones(size(X));
                    SF1_contribution(dist > upper_bound(K_index,1)) = 0;
                    SF1_contribution(dist < lower_bound(K_index,1)) = 0;
                    SF1 = SF1 + SF1_contribution;
                    SF2_contribution = K(K_index,5)*ones(size(X));
                    SF2_contribution(dist > upper_bound(K_index,1)) = 0;
                    SF2_contribution(dist < lower_bound(K_index,1)) = 0;
                    SF2 = SF2 + SF2_contribution;
                    
                    estimated_time_remaining = toc*((length(K(:,1))/K_index)-1)/60;
                    waitbar(K_index/length(K(:,1)),W,...
                    sprintf('Estimated Time Remaining: %3f minutes',estimated_time_remaining));
                end
            else
                for K_index = 1:length(K(:,1))
                    dist = sqrt((X-K(K_index,1)).^2+(Y-K(K_index,2)).^2+(Z-K(K_index,3)).^2);
                    V_contribution = K(K_index,4)*ones(size(X));
                    V_contribution(dist > upper_bound(K_index,1)) = 0;
                    V_contribution(dist < lower_bound(K_index,1)) = 0;
                    V = V + V_contribution;
                    estimated_time_remaining = toc*((length(K(:,1))/K_index)-1)/60;
                    waitbar(K_index/length(K(:,1)),W,...
                    sprintf('Estimated Time Remaining: %3f minutes',estimated_time_remaining));
                end
            end
            if complex
                V = SF1.^2 +SF2.^2;
            end
            this.screen_V = V;
            close(W);
            disp('Complete');
            disp(' ');
            fprintf('Total Calculation Time: %3f minutes\n',toc/60);
            dpi = resolution_width/6; %pixels/in ASSUMING 6 [in] screen width
            fprintf('Image dpi: %i \n',dpi);
            this.Plot_Kikuchi;
        end
        
        function K = Threshold_Data(this,threshold_magnitude,complex)
            if nargin < 2 || isempty(threshold_magnitude)
                threshold_magnitude = 3;
            end
            if nargin < 3 || isempty(complex)
                complex = this.complex;
            end
            K(:,1:3) = this.reciprocal_points;   %import and threshold the reciprocal points and intensity data
            if complex
                K(:,4:5) = this.RE_IM;
                logintensity = log10(K(:,4).^2 + K(:,5).^2);
                threshold = max(logintensity) - threshold_magnitude;
                K(logintensity<threshold,:) = [];
                magK = sqrt(K(:,1).^2 + K(:,2).^2 + K(:,3).^2);
                K(magK<0.2,:) = [];
            else
               if isempty(this.intensity)
                   this.Calculate_Intensity
               end
               % R(:,4:5) = this.intensity;
               
               K(:,4) = this.intensity(:,1);
               
               threshold = max(abs(K(:,4)))*10^(-threshold_magnitude);
               K(abs(K(:,4))<threshold,:) = [];
               % Filter out Direct Beam
               magK = sqrt(K(:,1).^2 + K(:,2).^2 + K(:,3).^2);
               K(magK<0.2,:) = [];
            end
        end
        
        function [X,Y,Z]= Screen_Vectors(this,resolutionwidth,zone,screendist)
            if nargin < 2 || isempty(resolutionwidth)
                resolutionwidth = 1024;
            end
            resolutionwidth = round(resolutionwidth);
            if nargin < 3 || isempty(zone)
                zone = [0,0,0];
            end
            if nargin < 4 || isempty(screendist)
                screendist = 3; %in
            end
            ilambda = 1/this.lambda;
            if zone(:) == 0
                p_inc = 2*pi/(resolutionwidth-1);
                a_inc = 2*pi/(resolutionwidth-1);
                ScreenX = zeros(resolutionwidth,resolutionwidth/2);
                ScreenY = zeros(resolutionwidth,resolutionwidth/2);
                ScreenZ = zeros(resolutionwidth,resolutionwidth/2);
                for p_index = 1:resolutionwidth
                    for a_index = 1:resolutionwidth/2
                        p_ang = (p_index-1)*p_inc;
                        a_ang = (a_index-1)*a_inc;
                        ScreenX(p_index,a_index) = ilambda*cos(p_ang)*sin(a_ang);
                        ScreenY(p_index,a_index) = ilambda*sin(p_ang)*sin(a_ang);
                        ScreenZ(p_index,a_index) = ilambda*cos(a_ang);
                    end
                end
                this.screen_X = ScreenX;
                this.screen_Y = ScreenY;
                this.screen_Z = ScreenZ;
                X = ScreenX;
                Y = ScreenY;
                Z = ScreenZ;
            else
                screenwidth=6; %in
                screenheight=4; %in
                resolutionheight=round(resolutionwidth*screenheight/screenwidth);
                [ScreenY,ScreenZ]=meshgrid(linspace(-screenwidth/2,screenwidth/2,resolutionwidth),...
                linspace(-screenheight/2,screenheight/2,resolutionheight));
                ScreenX = ones(size(ScreenY))*screendist;
                this.screen_X = ScreenX; this.screen_Y = ScreenY; this.screen_Z = ScreenZ;
                MAG = sqrt(ScreenX.^2+ScreenY.^2+ScreenZ.^2);
                X = ilambda*ScreenX./MAG;
                Y = ilambda*ScreenY./MAG;
                Z = ilambda*ScreenZ./MAG;
                screen_size = size(X);
                array_length = screen_size(1)*screen_size(2);
                screen_array(1,:)=reshape(X,[1,array_length]);
                screen_array(2,:)=reshape(Y,[1,array_length]);
                screen_array(3,:)=reshape(Z,[1,array_length]);
                zone_norm = sqrt(zone(1)^2 + zone(2)^2 + zone(3)^2);
                x = zone(1)/zone_norm;
                y = zone(2)/zone_norm;
                z = zone(3)/zone_norm;
                xy = sqrt(x^2 + y^2);
                a_rotator = [xy,0,-z;0,1,0;z,0,xy];
                p_rotator = [x/xy,-y/xy,0;y/xy,x/xy,0;0,0,1];
                screen_array = p_rotator*a_rotator*screen_array;
                X = reshape(screen_array(1,:),screen_size);
                Y = reshape(screen_array(2,:),screen_size);
                Z = reshape(screen_array(3,:),screen_size);
                % Normalize the Vectors and output them as X,Y,Z arrays, but
                % also save flat screen vectors to object properties for
                % display later
            end
        end
        
        function Plot_Kikuchi(this,X,Y,Z,V)
            if nargin < 2 || isempty(X)
                X = this.screen_X;
            end
            if nargin < 3 || isempty(Y)
                Y = this.screen_Y;
            end
            if nargin < 4 || isempty(Z)
                Z = this.screen_Z;
            end
            if nargin < 5 || isempty(V)
                V = this.screen_V;
            end
            
           if this.zone(:) == 0
               V = 100*V/max(max(V));
               figure1 = figure;
               axes1 = axes;
               surf(X,Y,Z,abs(V),'linestyle','none');
               set(figure1,'Color','black');
               campos([1,1,1]*10/this.lambda);
               % shading interp
               colormap(gray);
               grid(axes1,'off');
               set(axes1,'Color','black',...
                  'CLim',[0,50],...
                  'DataAspectRatio',[1,1,1],...
                  'XColor','black','YColor','black','ZColor','black',...
                  'Xlim',[-1/this.lambda,1/this.lambda],...
                  'Ylim',[-1/this.lambda,1/this.lambda],...
                  'Zlim',[-1/this.lambda,1/this.lambda]);
              axis vis3d
              zoom(1.75);
           else
               figure1 = figure;
               V = V/max(max(V));
               imshow(V);
               % set(gca,'Clim',[0,max(max]);
           end
        end
        
        function SAED_Calculation(this,zone,resolutionwidth,threshold_magnitude,sigma)
            if nargin < 2 || isempty(zone)
                zone = [1,0,0];
            end
            if nargin < 3 || isempty(resolutionwidth)
                resolutionwidth = 1024;
            end
            if nargin < 4 || isempty(threshold_magnitude)
                threshold_magnitude = 3;
            end
            if nargin < 5 || isempty(sigma)
                sigma = 0.5*this.reciprocal_spacing;
            end
            if this.complex
                this.Calculate_Intensity
            end
            ilambda = 1/this.lambda;
            this.zone = zone;
            [X,Y,Z] = this.Screen_Vectors(resolutionwidth,-(zone),50);
            zone_mag = sqrt(zone(1)^2 + zone(2)^2 + zone(3)^2);
            X = X+ilambda*zone(1)/zone_mag;
            Y = Y+ilambda*zone(2)/zone_mag;
            Z = Z+ilambda*zone(3)/zone_mag;
            R = this.Threshold_Data(threshold_magnitude);
            Ecenter = ilambda*zone/zone_mag;
            E_dist = sqrt((R(:,1)-Ecenter(1)).^2+(R(:,2)-Ecenter(2)).^2+(R(:,3)-Ecenter(3)).^2);
            R(abs(E_dist-ilambda)>sigma,:) = [];
            array_size = size(X);
            array_length = array_size(1)*array_size(2);
            X = reshape(X,array_length,1);
            Y = reshape(Y,array_length,1);
            Z = reshape(Z,array_length,1);
            V = zeros(array_length,1);
            tic
            W = waitbar(0,'Please Wait...','Name','SAED Calculation');
            for R_index = 1:length(R(:,1))
                dist = sqrt((X-R(R_index,1)).^2+(Y-R(R_index,2)).^2+(Z-R(R_index,3)).^2);
                V = V + normpdf(dist,0,sigma)*R(R_index,4);
                estimated_time_remaining = toc*((length(R(:,1))/R_index)-1)/60;
                waitbar(R_index/length(R(:,1)),W,...
                sprintf('Estimated Time Remaining: %1f minutes',estimated_time_remaining));
            end
%             Vmag = sqrt(X.^2 + Y.^2 + Z.^2);
%             V(Vmag < 0.5,:) = 0;
            V = reshape(V,array_size);
            this.SAED_Intensity = V;
            close(W);
            disp('Complete');
            disp(' ');
            fprintf('Total Calculation Time: %3f minutes\n',toc/60);
            this.Plot_SAED(zone,sigma,threshold_magnitude);
        end
        
        function Plot_SAED(this,zone,shell_thickness,magnitude_threshold)
            % method for creating crude SAED plots from diffraction data
            
            % shell_thickness is twice the maxiumum error bound for points
            % to be included in the Ewald sweep
            
            % threshold_magnitude is the number of orders of magnitude,
            % below which, intensity data will be excluded from the
            % calculation.
           if isempty(this.SAED_Intensity)
               if nargin < 2 || isempty(zone)
                   zone = [1,0,0]; % default to [1 0 0] zone
               end
               if nargin < 3 || isempty(shell_thickness) % Default shell thickness set to 1/2 the reciprocal sample lattice spacing
                   shell_thickness = this.reciprocal_spacing*0.5;
               end
               if nargin < 4 || isempty(magnitude_threshold) %How many orders of magnitude down from the maximum intensity do you want to threshhold your data at in this plot? Default is 2.
                   magnitude_threshold = 3;
               end
               if this.complex
                   this.Calculate_Intensity
               end
               Ecenter = zone/(sqrt(zone(1,1)^2 + zone(1,2)^2 + zone(1,3)^2)*this.lambda); % determine center of Ewald Sphere
               count = 1;
               V = zeros(1,4);
               for ii = 1:length(this.reciprocal_points(:,1)) % ignore all points that do not lie within shell thickness or are below intensity threshhold
                   if  abs(sqrt((Ecenter(1,1)-this.reciprocal_points(ii,1))^2 +...
                       (Ecenter(1,2)-this.reciprocal_points(ii,2))^2 +...
                       (Ecenter(1,3)-this.reciprocal_points(ii,3))^2)-(1/this.lambda)) < shell_thickness &&...
                       this.intensity(ii,2) > max(this.intensity(:,2))-magnitude_threshold
                      V(count,1:3) = this.reciprocal_points(ii,1:3);
                      V(count,4) = this.intensity(ii,2);
                      count = count+1;
                   end
               end
                figure1 = figure;
                axes1 = axes('Parent',figure1);
                set(gcf,'Color','black');
                scatter3(axes1,V(:,1),V(:,2),V(:,3),100,V(:,4),'.'); % Plot all non-ignored points in a scatter plot with various features
                colormap(gray);
                grid(axes1,'off');
                set(axes1,'CameraPosition',Ecenter,...
                    'Color','black',...
                    'CLim',[max(this.intensity(:,2))-magnitude_threshold,max(this.intensity(:,2))],...
                    'DataAspectRatio',[1,1,1],...
                    'XColor','black','YColor','black','ZColor','black',...
                    'Xlim',[min(this.reciprocal_points(:,1))*1.5,max(this.reciprocal_points(:,1))*1.5],...
                    'Ylim',[min(this.reciprocal_points(:,2))*1.5,max(this.reciprocal_points(:,2))*1.5],...
                    'Zlim',[min(this.reciprocal_points(:,3))*1.5,max(this.reciprocal_points(:,3))*1.5]);
           else
               figure;
               K = log10(this.SAED_Intensity);
               K_max = max(max(K));
               K_min = K_max - magnitude_threshold;
               imshow(K);
               set(gca,'Clim',[K_min,K_max]);
           end
        end
        
        function Calculate_Intensity(this)
            if this.complex && ~isempty(this.RE_IM)
               this.intensity(:,1) = (this.RE_IM(:,1).^2 + this.RE_IM(:,2).^2)/this.N_atoms; % calculate the intensity values from the structure factor
               this.intensity(:,2) = log10(this.intensity(:,1));
            else
                disp('The data is not complex or there is no structure factor data');
            end
        end
        
        function Kikuchi_Line_Scan(this,dpi,threshold_magnitude)
            
            if nargin < 2 || isempty(dpi)
                dpi = 4/(3*this.reciprocal_spacing*this.lambda);
            end
            if nargin < 3 || isempty(threshold_magnitude)
                threshold_magnitude = 3;
            end
            
            scan_zone = [1,0.5,0];
            scan_height = 4; %in
            screen_dist = 3; %in
            pixel_count = scan_height*dpi;
            zone_factor = screen_dist/sqrt(scan_zone(1)^2+scan_zone(2)^2+scan_zone(3)^2);
            
            X = zone_factor*ones(pixel_count,1)*scan_zone(1);
            Y = zone_factor*ones(pixel_count,1)*scan_zone(2);
            Z = (linspace(-scan_height/2,scan_height/2,pixel_count))';
            
            ilambda = 1/this.lambda;
            dK = this.reciprocal_spacing;
            pixel_norm = sqrt(X.^2+Y.^2+Z.^2);
            
            X = ilambda*X./pixel_norm;
            Y = ilambda*Y./pixel_norm;
            Z = ilambda*Z./pixel_norm;
            
            R = this.Threshold_Data(threshold_magnitude);
            V = zeros(size(X));
            
            R_mag = sqrt(R(:,1).^2 + R(:,2).^2 + R(:,3).^2);
            upper_bound = sqrt(ilambda^2 + (R_mag*dK/2));
            lower_bound = sqrt(ilambda^2 - (R_mag*dK/2));
            
            tic
            W = waitbar(0,'Please Wait...','Name','Kikuchi Line Scan Calculation');
           
            for R_index = 1:length(R(:,1))
                    dist = sqrt((X-R(R_index,1)).^2+(Y-R(R_index,2)).^2+(Z-R(R_index,3)).^2);
                    V_contribution = R(R_index,4)*ones(size(X));
                    V_contribution(dist > upper_bound(R_index,1)) = 0;
                    V_contribution(dist < lower_bound(R_index,1)) = 0;
                    V = V + V_contribution;
                    estimated_time_remaining = toc*((length(R(:,1))/R_index)-1)/60;
                    waitbar(R_index/length(R(:,1)),W,...
                    sprintf('Estimated Time Remaining: %3f minutes',estimated_time_remaining));
            end
            Z = Z.*pixel_norm/ilambda;
            close(W);
            disp('Complete');
            disp(' ');
            fprintf('Total Calculation Time: %3f minutes\n\n',toc/60);
            figure1 = figure;
            plot (Z,V);
            xlim([-0.1,0.1]);
            
        end
    end
    
end

