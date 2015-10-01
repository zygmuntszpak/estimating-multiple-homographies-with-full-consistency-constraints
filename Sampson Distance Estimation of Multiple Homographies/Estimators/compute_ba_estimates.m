%   Function: compute_ba_estimates
%
%   Given a series of corresponding points and several planes, this 
%   function computes homographies  so that they satisfy all possible
%   two-view consistency  constraints that arise due to the rigidity
%   of the scene. For example, two homographies corresponding to two 
%   facades of a building can be made fully consistent because the building
%   is rigid. However, if you somehow model the side views of two cars with
%   two different homographies, and these cars change their relative
%   positions then these homographies cannot be made fully consistent. 
%
%   Bundle adjustment is the gold-standard maximum likelihood estimation 
%   method. Suitable initial values are required for the latent variables 
%   that the bundle adjustment operates on, otherwise the method
%   will converge to a poor solution. Good initial values can be found
%   using an iterative 'compute_aml_estimates'  scheme or by using the 
%   so-called direct Chojnacki initialisation method. 
%   The Chojnacki initialisation method is much faster. 
%
%   Parameters:
%
%      sceneData 	   			   - a struct containing noisy 
%                                    corresponding points between two views
%
%      listOfInitialH  			   - initial homography estimates, either 
%                                    DLT or FNS based
%
%      useChojnackiInitialisation  - simple flag variable, toggling between
%                                    using a Chojnacki initialisation to 
%                                    initialise the required latent 
%                                    variables or Chen's initialisation 
%                                    method. Chojnacki's method works
%									 when there are two or more 
%                                    homographies in the scene, while
%									 Chen's requires at least three
%                                    homographies.
% 
%										
%
%   Returns: A cell array of fully consistent homographies.
% 
%
%   See Also:
%
%  compute_aml_estimates
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

function  [listOfEstimatedH, diagnostic] = compute_ba_estimates( ...
                    sceneData, listOfInitialH, useChojnackiInitialisation )
    
    % data structure used to keep track of diagnostic information relating
    % to the performance of the optimization method (e.g. running time etc)
    diagnostic.iterations = 0;
    diagnostic.residualHistory = [];
    diagnostic.runningTime = 0;

    % transfer the pts in all views and the list of initial homographies to
    % a common coordinate system for numerical stability
    [sceneDataJointlyNormalised, listOfInitialH, T, Tp] = ...
                      normalise_joint_transform(sceneData,listOfInitialH);
    
    
    numOfH = length(sceneData);
    listOfEstimatedH = cell(1,numOfH);

    if (useChojnackiInitialisation)
      eta = initialise_eta_chojnacki(listOfInitialH);
    else
      eta = initialise_eta_peichen_true(listOfInitialH);  
    end
    
	numberOfPrincipalParameters = length(eta);
	
    % Append points in first view for all planes to the latent variable 
    % vector eta
	% recall that bundle adjustment needs these 'nuisance parameters'
    for i=1:numOfH
        pts = sceneDataJointlyNormalised(i).ptsInViewOne;
        eta = [eta(:) ; pts(1,:)' ; pts(2,:)'];
    end
    
    %'TolFun',1e-8,'TolX',1e-8
     options = optimset('FinDiffType ','central','Diagnostics','on',...
         'Display','iter','Algorithm',{'levenberg-marquardt',10},...
                                            'TolFun',1e-10,'TolX',1e-10,...
                              'Jacobian','off','DerivativeCheck','off',...
                              'OutputFcn', @myoutput,...
                              'MaxIter',100000);
          
    f = @(bestEta)cost_function_bundleadjustment(...
                                       bestEta,sceneDataJointlyNormalised);
    
    % keeps track of the residuals for each iterations of the optimisation
    residualHistory = [];
    tic % start timer                               
    [etaFinal,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                            lsqnonlin(f,eta,[],[],options);
    runningTime = toc; % end timer                                  
    % Nested function whose sole purpose is to keep track of the
    % AML cost function at each iteration. The values will be used later
    % to analysie the convergence properties of our algorithm
    function stop = myoutput(~,optimvalues,state)
        stop = false;
        if state == 'iter'
          residualHistory = [residualHistory; optimvalues.resnorm];
        end
    end

    % etaFinal contains the final values of the latent parameters encoding 
    % the consistent homographies and all the points in the first view.
    % We want to extract the homographies from eta
	etaFinal = etaFinal(1:numberOfPrincipalParameters);
    listOfEstimatedH = convertEtaToHomographyCellArray(etaFinal);
    
    for i = 1:length(listOfEstimatedH)
        listOfEstimatedH{i} =  Tp \ listOfEstimatedH{i} * T;
    end
    
    % store the diagnostic information 
    diagnostic.iterations = output.iterations;
    diagnostic.residualHistory = residualHistory;
    diagnostic.runningTime = runningTime;

end

%   Function: convertEtaToHomographyCellArray
%
%   Converts a vector of latent variables encoding a collection of fully 
%   consistent homographies  into individual homography matrices.
%
%   The structure of the eta parameter is slightly different from the
%   notation used in the paper. 
%   
%   The eta vector is : [a, b, v_1,w_1, ... , v_n, w_n], whereas in the
%   paper we structure it as [w, b, v_1, ... , v_n, w_1, ...,w_n]; 
%   Keep in mind that this has implications for the overall structure
%   of the Jacobian matrix.
%
%   Parameters:
%
%      eta 	   			     - vector of latent variables encoding 
%                              the consistent homographies
%
%										
%
%   Returns: A cell array containing fully consistent homography matrices
% 
function listOfOptimisedValidHomographies = ...
                                      convertEtaToHomographyCellArray(eta)

% there are 4*I+12 parameters for I planes
numberOfPlanes = (length(eta) - 12) / 4;
listOfOptimisedValidHomographies = cell(1,numberOfPlanes);

% extract first latent variable
a = eta(1:9);
% convert length-9 vector into 3x3 matrix
A = reshape(a,3,3);
% extact second latent variable
b = eta(10:12);
index = 13;
% Extract remaining latent variables
for i = 1:numberOfPlanes
    vi = eta(index:index+2);
    wi = eta(index+3);
    index = index + 4;
    % construct the i'th homography using the latent variables
    listOfOptimisedValidHomographies{i} = (wi * A + b * vi'); 
end

end

%   Function: cost_function_bundleadjustment
%
%   Computes the quality of the homographies encoded by the latent 
%   variables in the vector eta. The errors are the 'reprojection' errors
%   for each point correspondence. 
%
%   Parameters:
%
%      eta 	     - vector of latent variables encoding the consistent 
%                  homographies
%     
%      sceneData - A struct containing noisy point correspondences
%                  between two views
%
%										
%
%   Returns: A vector of reprojection errors
% 
function errors = cost_function_bundleadjustment(eta, sceneData )
    numberOfPlanes = length(sceneData);
    listOfEstimatedH = cell(1,numberOfPlanes);
    listOfEstimatedP = cell(1,numberOfPlanes);

    a = eta(1:9);
    % vertical concatenating the rows to form 3 into 3
    % rotation matrix
    A = reshape(a,3,3);
    % extract the translation of the camera
    b = eta(10:12);

    index = 13;
    % Extract the plane normals and distance from origin
    for i = 1:numberOfPlanes
        vi = eta(index:index+2);
        wi = eta(index+3);
        index = index + 4;
        listOfEstimatedH{i} = (wi * A + b * vi');  

    end
    
    % Extract the estimated points
    indexWherePointsStart = index;
    totalNumberOfPoints = 0;
    for k = 1:numberOfPlanes
        numberOfPoints = length(sceneData(k).ptsInViewOne);
        estimatedPointsInPlaneK = ...
            eta(indexWherePointsStart:indexWherePointsStart+...
                                                    (numberOfPoints*2)-1);
        modelCorrespondingPointsMatrix = reshape(...
                              estimatedPointsInPlaneK,numberOfPoints,2)';
        listOfEstimatedP{k} = modelCorrespondingPointsMatrix;
        %listOfEstimatedP{k} - sceneData(k).ptsInViewOne;
        indexWherePointsStart = indexWherePointsStart+(numberOfPoints*2);
        totalNumberOfPoints = totalNumberOfPoints + (numberOfPoints*2);
        
    end
    %errors = zeros(1,totalNumberOfPoints);
    errors = [];
    % Compute reprojection error
    for k = 1:numberOfPlanes
        H = listOfEstimatedH{k};
        x = listOfEstimatedP{k};
        n = length(x);
        % Generate the corresponding points xp in image 2
        xp = H * [ x; ones( 1, n ) ];
        xp = xp ./ repmat( xp(3,:), 3, 1 );
        xp = xp(1:2,:);
        1;
        
        errorsViewOne = sceneData(k).ptsInViewOne-x;
        errorsViewTwo = sceneData(k).ptsInViewTwo-xp; 
        combinedErrors = [errorsViewOne; errorsViewTwo];
        
        errors = [errors(:) ; combinedErrors(:) ];
        1;

    end
    
end

function [sceneData, listOfInitialH, T, Tp] = ...
                       normalise_joint_transform(sceneData, listOfInitialH)
    numOfH = length(sceneData);

    
    % collect all the pts in view one into a matrix
    % collect all the pts in view two into a matrix
    x = sceneData(1).ptsInViewOne;
    xp = sceneData(1).ptsInViewTwo;
    for i = 2:numOfH
      x = horzcat(x,sceneData(i).ptsInViewOne);
      xp = horzcat(xp,sceneData(i).ptsInViewTwo);
    end
    
    nn = length(x);
    x = [ x; ones( 1, nn ) ];
    xp = [ xp; ones( 1, nn ) ];
    T = normalize_transform( x );
    Tp = normalize_transform( xp );
    
    % convert the points and homographys in each plane into the new
    % joint coordinate system    
    for i = 1:numOfH
        %listOfInitialH{i} = Tp \ listOfInitialH{i} * T;
        listOfInitialH{i} = Tp * listOfInitialH{i} / T;
        x = sceneData(i).ptsInViewOne;
        xp = sceneData(i).ptsInViewTwo;
        n = length(x);
        x = [ x; ones( 1, n ) ];
        xp = [ xp; ones( 1, n ) ];
        
        xx =  T * x;
        xpp = Tp * xp;
        
        xx = xx ./ repmat( xx(3,:), 3, 1 );
        xpp = xpp ./ repmat( xpp(3,:), 3, 1 );
        sceneData(i).ptsInViewOne = xx(1:2,:);
        sceneData(i).ptsInViewTwo = xpp(1:2,:);       
    end

end

function T = normalize_transform( x )

  % Transform taking x's centroid to the origin

  Ttrans = [ 1 0 -mean( x(1,:) ) ; 0 1 -mean( x(2,:) ) ; 0 0 1 ];

  % Calculate appropriate scaling factor

  x = Ttrans * x;
  lengths = sqrt( sum( x(1:2,:).^2 ));
  s = sqrt(2) / mean(lengths);

  % Transform scaling x to an average length of sqrt(2)

  Tscale = [ s 0 0 ; 0 s 0 ; 0 0 1 ];

  % Compose the transforms

  T = Tscale * Ttrans;
end

