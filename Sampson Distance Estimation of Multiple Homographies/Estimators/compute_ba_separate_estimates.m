%   Function: compute_ba_separate_estimates
%
%   Given a series of corresponding points and several planes, this 
%   function computes homographies  for each of the planes by minimising
%   the 'reprojection' error, but without enforcing inter homography 
%   constraints. In other words, it applies what would be considered the
%   maximum likelihood estimate for a single homography matrix. However, 
%   because there are multiple planes in the scene it is no longer truly
%   a maximum likelihood procedure because it does not enforce inter 
%   homography constraints.
%
%
%   Parameters:
%
%      sceneData 	   			   - a struct containing noisy 
%                                    corresponding points between two views
%
%      listOfInitialH  			   - initial homography estimates, either 
%                                    DLT or FNS based
%
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

function  [listOfEstimatedH, diagnostic] = ...
                  compute_ba_separate_estimates(sceneData, listOfInitialH)
    
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
    
    eta = [];
    % vectorise homography matrices and pack them into one long parameter
    % vector
    for i = 1:numOfH
        tempH = listOfInitialH{i};
        eta = [eta(:); tempH(:)];
    end
    
	numberOfPrincipalParameters = length(eta);
	
    % Append points in first view for all planes to the parameter 
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

    % etaFinal contains the final values of the homography parameters 
    % We want to extract the homographies from eta
	etaFinal = etaFinal(1:numberOfPrincipalParameters);
    k = 1;
    for i = 1:9:numOfH*9
        tempH = etaFinal(i:i+8);
        listOfEstimatedH{k} = reshape(tempH,3,3); 
        k = k +1;
    end

    
    for i = 1:length(listOfEstimatedH)
        listOfEstimatedH{i} =  Tp \ listOfEstimatedH{i} * T;
    end
    
    % store the diagnostic information 
    diagnostic.iterations = output.iterations;
    diagnostic.residualHistory = residualHistory;
    diagnostic.runningTime = runningTime;

end



%   Function: cost_function_bundleadjustment
%
%   Computes the quality of the homographies encoded by the latent
%   variables in the vector eta.  The errors are the 'reprojection' errors
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

    i = 1;
    for k = 1:9:numberOfPlanes*9
        tempH = eta(k:k+8);
        listOfEstimatedH{i} = reshape(tempH,3,3); 
        i = i + 1;
    end    

    index = (numberOfPlanes*9)+1;   
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

