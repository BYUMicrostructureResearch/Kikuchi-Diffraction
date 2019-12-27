clc; clearvars; close all
% Version 1.3.0 - 22 Dec 2015 - Last edits by Shawn P. Coleman
% Version 1.2.5 - 27 Oct 2015 - Last edits by Adam D. Herron

%% User input

dpi = 256;
screendist = 3;   %in
screenwidth = 6;  %in
screenheight = 4; %in
beam_radius=0.2;  
threshold_magnitude = 3.5;  %% Note, there is also threshold option on the fix

filename='test_Fe.0.xyz';
zone=[1 1 1];



%% Import intensity data from fix_saed_xyz

fid=fopen(filename);
lambda=textscan(fid,'%*s %f',1);
lambda=cell2mat(lambda);
reciprocal_spacing=textscan(fid,'%*s %f %f %f',1);
reciprocal_spacing=min(cell2mat(reciprocal_spacing));
data=textscan(fid,'%f %f %f %f');
fclose(fid);
data=cell2mat(data);

%% Screen Vectors
resolution_width = round(dpi*screenwidth);
ilambda = 1/lambda;

if zone(:) == 0 % Sphere
  p_inc = 2*pi/(resolution_width-1);
  a_inc = 2*pi/(resolution_width-1);
  ScreenX = zeros(resolution_width,resolution_width/2);
  ScreenY = zeros(resolution_width,resolution_width/2);
  ScreenZ = zeros(resolution_width,resolution_width/2);
  for p_index = 1:resolution_width
    for a_index = 1:resolution_width/2
      p_ang = (p_index-1)*p_inc;
      a_ang = (a_index-1)*a_inc;
      ScreenX(p_index,a_index) = ilambda*cos(p_ang)*sin(a_ang);
      ScreenY(p_index,a_index) = ilambda*sin(p_ang)*sin(a_ang);
      ScreenZ(p_index,a_index) = ilambda*cos(a_ang);
    end
  end
  X = ScreenX;
  Y = ScreenY;
  Z = ScreenZ;
else
  resolution_height=round(dpi*screenheight);
  [ScreenY,ScreenZ]=meshgrid(linspace(-screenwidth/2,screenwidth/2,resolution_width),linspace(-screenheight/2,screenheight/2,resolution_height));
  ScreenX = ones(size(ScreenY))*screendist;

  MAG = sqrt(ScreenX.^2+ScreenY.^2+ScreenZ.^2);
  X = ilambda*ScreenX./MAG;
  Y = ilambda*ScreenY./MAG;
  Z = ilambda*ScreenZ./MAG;
  screen_size = size(X);
  array_length = resolution_width*resolution_height;
  screen_array(1,:)=reshape(X,1,[]);
  screen_array(2,:)=reshape(Y,1,[]);
  screen_array(3,:)=reshape(Z,1,[]);
  zone_norm = norm(zone);
  x = zone(1)/zone_norm;
  y = zone(2)/zone_norm;
  z = zone(3)/zone_norm;

  xy = sqrt(x^2 + y^2);  
  a_rotator = [xy,0,-z;0,1,0;z,0,xy];
  if (xy~=0) 
    p_rotator = [x/xy,-y/xy,0;y/xy,x/xy,0;0,0,1];
  else  % resolve the issue when x=y=xy=0
    x_xy=sqrt(2)/2; y_xy=sqrt(2)/2;
    p_rotator = [x_xy,-y_xy,0;y_xy,x_xy,0;0,0,1];
  end
  
  screen_array = p_rotator*a_rotator*screen_array;
  X = reshape(screen_array(1,:),screen_size);
  Y = reshape(screen_array(2,:),screen_size);
  Z = reshape(screen_array(3,:),screen_size);
  % Normalize the Vectors and output them as X,Y,Z arrays, but
  % also save flat screen vectors to object properties for
  % display later
end

%% Threshold & Filter Beam

K = data;   %import and threshold the reciprocal points and intensity data

threshold = max(abs(K(:,4)))*10^(-threshold_magnitude);
K(abs(K(:,4))<threshold,:) = [];

K_mag = sum(K(:,1:3).^2')';
K(K_mag<beam_radius^2,:) = [];
K_mag(K_mag<beam_radius^2,:) = [];

%% Collection size
% dK = reciprocal_spacing; % Actual dK
% dK = 4*ilambda*((1+((R_mag.^2)*(lambda^2)/4)).^(-3/2))/(screendist*dpi); %Alternate dK_eff
dK = 4*ilambda/sqrt(((dpi^2)*(screendist^2)+1)); % dK_eff (lower limit on dK for specified dpi)

if dK < reciprocal_spacing
  dK = reciprocal_spacing; % minimum bandwidth for a single sample point is 1-2 pixels
end
            
if zone(:) == 0
  dK = 2*dK;
end

% upper_bound = sqrt(ilambda^2 + (R_mag*dK/2)); % Actual dK
% lower_bound = sqrt(ilambda^2 - (R_mag*dK/2)); % Actual dK
% upper_bound = sqrt((ilambda^2)-K_mag); % All the way to theta == 0
upper_bound2 = ((ilambda^2)+K_mag.*dK/2); % Req dK
lower_bound2 = ((ilambda^2)-K_mag.*dK/2); % Req dK

%% Kikuchi Calculation 
tic
step=0;
W = waitbar(0,'Please Wait...','Name','Classic Kikuchi Calculation');
V = zeros(size(X));
for K_index = 1:length(K(:,1))
  dist2 = ((X-K(K_index,1)).^2+(Y-K(K_index,2)).^2+(Z-K(K_index,3)).^2);                 % expensive
  V= V+[dist2 < upper_bound2(K_index,1) & dist2 > lower_bound2(K_index,1)]*K(K_index,4);        % expensive
 
  if K_index/length(K(:,1))>step
    estimated_time_remaining = toc*((length(K(:,1))/K_index)-1)/60;
    waitbar(K_index/length(K(:,1)),W,...
      sprintf('Estimated Time Remaining: %3.3f minutes',estimated_time_remaining));
    step=step+0.01;
  end
  
end
close(W);

disp('Complete');
disp(' ');
fprintf('Total Calculation Time: %3f minutes\n',toc/60);
dpi = resolution_width/6; %pixels/in ASSUMING 6 [in] screen width
fprintf('Image dpi: %i \n',dpi);


%% Plot Kikuchi 
if zone(:) == 0
  V = 100*V/max(max(V));
  figure1 = figure;
  axes1 = axes;
  surf(X,Y,Z,abs(V),'linestyle','none');
  set(figure1,'Color','black');
  campos([1,1,1]*10/lambda);
  % shading interp
  colormap(gray);
  grid(axes1,'off');
  set(axes1,'Color','black',...
    'CLim',[0,50],...
    'DataAspectRatio',[1,1,1],...
    'XColor','black','YColor','black','ZColor','black',...
    'Xlim',[-1/lambda,1/lambda],...
    'Ylim',[-1/lambda,1/lambda],...
    'Zlim',[-1/lambda,1/lambda]);
  axis vis3d
  zoom(1.75);
else
  V = V/max(max(V));
  figure1 = figure;
  imshow(V);
% % (note, exportfig is a function from matlab file exchange)  
%   expfig=['export_fig ' filename(1:end-6) '_' num2str(zone(1)) '_' num2str(zone(2)) '_' num2str(zone(3))  ' -png'];
%   eval(expfig)
%   close all

  figure2 = figure;
  imshow(abs(V-0.5));
%   expfig=['export_fig ' filename(1:end-6) '_' num2str(zone(1)) '_' num2str(zone(2)) '_' num2str(zone(3))  '_gray -png'];
%   eval(expfig)
%   close all
end

% save([filename(1:end-6) '_' num2str(zone(1)) '_' num2str(zone(2)) '_' num2str(zone(3)) '.mat'])

