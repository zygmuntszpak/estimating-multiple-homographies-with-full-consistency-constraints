%   Function: generate_groundtruth_scene_new
%
%   Generates a random planar scene consisting of a specified number of 
%   planes. Additionally, two cameras facing the planes are generated and 
%   a specified number of 3D points on the 3D planes are projected onto two
%   images using  Then the groundt truth homographies that map 
%   corresponding points between the two views are computed.
%
%   Parameters:
%
%      numH               - number of 3D planes, and hence the number of
%							homographies to generate
%
%      nPoints            - number of corresponding points 
%
%
%
%   Returns:
%
%     A sceneData struct that contains the homographies and corresponding
%     points from view one to view two.
%
%   See Also:
%
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 16/9/2014

function sceneData = generate_groundtruth_scene_new(numOfH,nPoints)
format long g

% data structures for homographies, points in the first view and points
% in the second view
H_list = cell(1,numOfH);
x_list = cell(1,numOfH);
xp_list = cell(1,numOfH);

% we will keep generating random scenes until all homographies have
% at least half the requested number of corresponding points
minPoints = nPoints*0.5;
        
for i = 1:numOfH
        
    % a homography is considered reasonable if the number of visible
    % points in the second view is >= minN
    isReasonableHomography = 0;
    
    while (~isReasonableHomography)
        
        % focal lenght
        f = 250;
        
        % intersection of optical axis and image
        px = 0;
        py = 0;
        
        % calibration matrix for first camera
        K1 = [f,0,px;
            0, f, py ;
            0, 0, 1 ];
        
        
        % calibration matrix for second camera
        K2 = [f+100,0,px;
            0, f+100, py;
            0, 0, 1];
        
        
        % rotation matrix of the first camera
        alpha = 10 * (pi/180);
        beta = 5 * (pi/180);
        R1 = [1,0,0;
            0, cos(alpha), -sin(alpha);
            0, sin(alpha), cos(alpha)] *  ...
            [cos(beta), 0, sin(beta);
            0,         1,        0 ;
            -sin(beta), 0 , cos(beta)];
        
        % rotation matrix of the second camera
        alpha2 = -5 * (pi / 180);
        beta2 = 10 * (pi / 180);
        R2 = [1,0,0 ;
            0, cos(alpha2), -sin(alpha2);
            0, sin(alpha2), cos(alpha2)] * ...
            [cos(beta2), 0 , sin(beta2) ;
            0,            1,          0;
            -sin(beta2), 0 , cos(beta2)];
        
        % translation matrix for the first camera
        tCamera1 = [-150; 0 ;-100];
        
        % translation matrix for the second camera
        tCamera2 = [200; 100; -100];
                
        
        % unit normal originaly pointing toward camera
        normal = cross([10 10 0], [15 5 0]);
        unitNormal = (normal/norm(normal))';
        
        
        % select a random distance for this plane as measured
        % from the origin
        initialDistanceOfPlaneFromOrigin = randi(650) + 650;
        p1 = [0 0 initialDistanceOfPlaneFromOrigin];
        
        % randomly rotate the normal so that it no-longer directly faces
        % the camera
        if (rand(1,1) >= 0.5)
            alpha1 =  (rand(1,1) * 30)  * (pi/180);
        else
            alpha1 =  (rand(1,1) * -30)  * (pi/180);
        end
        
        if (rand(1,1) >= 0.5)
            beta1 =  (rand(1,1) * 30)  * (pi/180);
        else
            beta1 =  (rand(1,1) * -30)  * (pi/180);
        end
        
        R = get_rotation_matrix(alpha1,beta1);
        unitNormal = R * unitNormal;
        unitNormal = unitNormal / norm(unitNormal);
        
        % make sure that we are consitent with our choice of what we 
        % consider "further away" from the origin (e.g. positive values
        % further away)
        distanceOfPlaneFromOrigin = -1*dot(p1,unitNormal);
        
        
        % Latent variables for the homography matrices...
        A = K2*R2*inv(R1)*inv(K1);
        b = K2*R2*(tCamera1-tCamera2);
        w = unitNormal'*tCamera1 - distanceOfPlaneFromOrigin;
        v = inv(K1)'*inv(R1)'*unitNormal;
        
        % Construct a geometrically consistent homography for the scene
        H = w*A + b*v';
        
        % we fix our image canvas to be of dimension 500 x 500 pixels
        % running from -250 to 250, with the origina at the centre of
        % the canvas. 
        boundaryEndC = 250;
        boundaryStartC = -250;
        boundaryEndR = 250;
        boundaryStartR = -250;
        
        % we would like the corresponding points associated with 
        % homographies to span spatial regions of varying sizes.
        % to achieve this, we first sample points inside random
        % sized rectangle in the centre of the image canvas
        regionBoundaryEndC =  randi([40,120],1);
        regionBoundaryStartC = -regionBoundaryEndC;
        regionBoundaryEndR = randi([40,120],1);
        regionBoundaryStartR = -regionBoundaryEndR;
         
        % distance to border of image
        maxShift = 250 - regionBoundaryEndC;        
        
        % we then wish to randomly shift the points 
        randomShiftC = randi([-maxShift,maxShift],1);
        randomShiftR = randi([-maxShift,maxShift],1);
   
        % assuming a square aspect ratio of image, the desired sampling
        % can be achieved with the following line of code
        x = randi([regionBoundaryStartC,regionBoundaryEndC],[2,nPoints])...
            + repmat([randomShiftC,randomShiftR]',1,nPoints);
        
        % now we transfer the points to the second view using the
        % homography
        xp = H * [ x; ones( 1, nPoints ) ];
        xp = xp ./ repmat( xp(3,:), 3, 1 );
        
        % finally we count how many of the transfered points are still
        % inside the bounds of the image
        
        ind = find(xp(1,:) <= boundaryEndC);
        % drop the homogenous coordinate and keep points that survive the 
        % first condition
        xp = xp(1:2,ind);
        x = x(1:2,ind);
        ind = find(xp(1,:) >= boundaryStartC);
        % keep points that survive the second condition
        xp = xp(1:2,ind);
        x = x(1:2,ind);
        ind = find(xp(2,:) <= boundaryEndR);
        xp = xp(1:2,ind);
        x = x(1:2,ind);
        ind = find(xp(2,:) >= boundaryStartR);
        
        %fprintf('number of points found: %f \n',length(ind))
        if (length(ind) >= minPoints)
            isReasonableHomography = 1;
            x = x(1:2,ind);
            xp = xp(1:2,ind);            
        end
    end
    
    H_list{i} = H;
    x_list{i} = x;
    xp_list{i} = xp;
    sceneData = struct('homographies',H_list,'ptsInViewOne',x_list,...
                                               'ptsInViewTwo',xp_list);
end


end

function R = get_rotation_matrix(alpha,beta)
% All rotations are counterclockwise about the respective axes and relative
% to the camera coordinate system.
% For each view, we first apply a rotation about the y-axis followed
% by one about the x-axis (there is no z-axis rotation).
% By composition of the two transformation we have the rotation matrix for
% the left camera
R = [1 0 0 ; 0 cos(alpha) -sin(alpha) ;
    0 sin(alpha) cos(alpha)] * [cos(beta) 0 sin(beta) ;
    0 1 0; -sin(beta) 0 cos(beta)];
end



