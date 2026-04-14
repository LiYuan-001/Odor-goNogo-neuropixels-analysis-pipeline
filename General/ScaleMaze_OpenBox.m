function [x_rotated_scaled,y_rotated_scaled, x_center,y_center,maze,scaling_for_later] = ScaleMaze_OpenBox(x, y, coord, target_maze,scaleStat,scale) 

ind = isnan(x);
if scaleStat == 0
    % % Scale and rotate maze if required
    % %-----------------------------------------------------
    % % Get Angle; Set center
    % %-----------------------------------------------------
    % 
    % theta = coord.theta_rad;
    % 
    % % choose a point which will be the center of rotation
    % x_center = mean(coord.line(:,1));
    % y_center = mean(coord.line(:,2));
    % 
    % % create a matrix which will be used later in calculations
    % center = repmat([x_center; y_center], 1, length(x));
    % 
    % %angle = %rotation_matrix = repmat([
    % % define a 60 degree counter-clockwise rotation matrix %theta = pi/3;       % pi/3 radians = 60 degrees
    % 
    % R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    % 
    % %-----------------------------------------------------
    % % Rotate running path
    % %-----------------------------------------------------
    % v = [x;y];          % create a matrix to apply
    % s = v - center;     % shift points in the plane so that the center of rotation is at the origin
    % so = R*s;           % apply the rotation about the origin
    % x_rotated = so(1,:);
    % y_rotated = so(2,:);
    % 
    % % Taking this line out will let your maze be centered
    % % vo = so + center;   % shift again so the origin goes back to the desired center of rotation                     
    % % this can be done in one line as:
    % % vo = R*(v - center) + center
    % 
    % %-----------------------------------------------------
    % % Now rotate poly
    % %-----------------------------------------------------
    % 
    % NE = [coord.NE(1); coord.NE(2)];NW = [coord.NW(1); coord.NW(2)];
    % SE = [coord.SE(1); coord.SE(2)];SW = [coord.SW(1); coord.SW(2)];
    % maze.NE_rot = R*(NE - [x_center; y_center]); %+ [x_center; y_center] %put either back to origin, or taking second part out, new center is 0
    % maze.NW_rot = R*(NW - [x_center; y_center]); %+ [x_center; y_center]
    % maze.SE_rot = R*(SE - [x_center; y_center]); %+ [x_center; y_center]
    % maze.SW_rot = R*(SW - [x_center; y_center]); %+ [x_center; y_center]
    % maze.xCenter = 0; maze.yCenter = 0;
    % % pick out the vectors of rotated x- and y-data

    %-----------------------------------------------------
    % Scaling
    %-----------------------------------------------------
    target(1) = target_maze(1); 
    target(2) = target_maze(2); 


    % 1: Scaling Arms proportionally in x and in y along the entire axis
    positive_x = find(x_rotated>0);
    negative_x = find(x_rotated<0);
    x_scaling(1) =(target(1) / (mean([maze.NE_rot(1) maze.SE_rot(1)])));
    x_rotated_scaled(positive_x) = x_rotated(positive_x) .* x_scaling(1); 
    x_scaling(2) =(target(1) / (mean([maze.NW_rot(1) maze.SW_rot(1)])));
    x_rotated_scaled(negative_x) = -x_rotated(negative_x) .* x_scaling(2);
    

    y_scaling = (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)]));
    y_rotated_scaled = y_rotated .* y_scaling;

    % 2: thes eare the dots from the GUI just for orientation
    maze.NE_rot_scaled(1) = maze.NE_rot(1) .* target(1) / (mean([maze.NE_rot(1) maze.SE_rot(1)]));maze.SE_rot_scaled(1) = maze.SE_rot(1) .* target(1) / (mean([maze.NE_rot(1) maze.SE_rot(1)])); 
    maze.NW_rot_scaled(1) = -maze.NW_rot(1) .* target(1) / (mean([maze.NW_rot(1) maze.SW_rot(1)]));maze.SW_rot_scaled(1) = -maze.SW_rot(1) .* target(1) / (mean([maze.NW_rot(1) maze.SW_rot(1)]));
    maze.NE_rot_scaled(2) = maze.NE_rot(2)* (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)])); maze.NW_rot_scaled(2) = maze.NW_rot(2) * (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)]));
    maze.SE_rot_scaled(2) = maze.SE_rot(2)* (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)])); maze.SW_rot_scaled(2) = maze.SW_rot(2) * (target(2) / mean([maze.NW_rot(2)-maze.SW_rot(2); maze.NE_rot(2)-maze.SE_rot(2)]));
    
    plot(x_rotated_scaled,y_rotated_scaled, 0, 0, 'ko')
    hold on
    scatter([maze.SW_rot_scaled(1), maze.SE_rot_scaled(1), maze.NW_rot_scaled(1), ...
        maze.NE_rot_scaled(1)], [maze.SW_rot_scaled(2), maze.SE_rot_scaled(2), maze.NW_rot_scaled(2), maze.NE_rot_scaled(2)], 'g')
    %-----------------------------------------------------
    % Scaling for Sleep and open box
    %-----------------------------------------------------
    scaling_for_later.xaxis = mean(abs(x_scaling));
    scaling_for_later.yaxis = mean(abs(y_scaling));
    
else % while it is open box and sleep box
    scaling_for_later = scale;
    
    x_center = (min(x)+max(x))/2;
    y_center = (min(y)+max(y))/2;
%     x = (x - x_center);
%     y = (y - y_center);
    
    % Scale and rotate maze if required
    %-----------------------------------------------------
    % Get Angle; Set center
    %-----------------------------------------------------

    theta = coord.theta_rad;

    % create a matrix which will be used later in calculations
    center = repmat([x_center; y_center], 1, length(x));

    %angle = %rotation_matrix = repmat([
    % define a 60 degree counter-clockwise rotation matrix %theta = pi/3;       % pi/3 radians = 60 degrees

    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    %-----------------------------------------------------
    % Rotate running path
    %-----------------------------------------------------
    v = [x;y];          % create a matrix to apply
    s = v - center;     % shift points in the plane so that the center of rotation is at the origin
    so = R*s;           % apply the rotation about the origin
    x_rotated = so(1,:);
    y_rotated = so(2,:);
    
    
    x_rotated_scaled = x_rotated * scale.xaxis;
    y_rotated_scaled = y_rotated * scale.yaxis;
    
    plot(x_rotated_scaled,y_rotated_scaled);axis equal
    hold on
    maze.center = [mean(x_rotated_scaled) mean(y_rotated_scaled)];
end

x_rotated_scaled(ind) = NaN;
y_rotated_scaled(ind) = NaN;
