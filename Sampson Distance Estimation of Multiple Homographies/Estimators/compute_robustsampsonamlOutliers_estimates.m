%   Function: compute_robustsampsonaml_estimates
%
%   Given a collection of homographies estimated using DLT or FNS together
%   with corresponding points this function 'upgrades' the homographies so
%   that they satisfy all possible two-view consistency constraints that
%   arise due to the rigidity of the scene. For example, two homographies 
%   corresponding to two facades of a building can be made fully consistent
%   because the building is rigid. However, if for some reason you model 
%   the side views of two cars with two different homographies, and these
%   cars change their relative positions between two views, 
%   then these homographies cannot be made fully consistent.
%   under the assumed model.
%
%   Parameters:
%
%      sceneData 	   			   - a struct containing noisy 
%                                    corresponding points between two views
%
%      listOfInitialH  			   - initial homography estimates, 
%                                    either DLT or FNS based
%
%      useChojnackiInitialisation  - simple flag variable, toggling between
%                                    using a Chojnacki nitialisation to 
%                                    initialise the required latent
%                                    variables or Chen's initialisation
%                                    method. Chojnacki's method works
%									 when there are two or more 
%                                    homographies in the scene, while
%									 Chen's requires at least three 
%                                    homographies.
%
%      useFNSCovariances           - simple flag variable, toggling between
%                                    using FNS covariances
%								     or DLT covariances.   
%										
%
%   Returns: A cell array of fully consistent homographies.
% 
%
%   See Also:
%
%  compute_ba_estimates
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

function [listOfEstimatedH, diagnostic] = ...
                            compute_robustsampsonamlOutliers_estimates(...
                sceneData, listOfInitialH, useChojnackiInitialisation,...
                                                        useFNSCovariances,...
                                                        threshold)

    % data structure used to keep track of diagnostic information relating
    % to the performance of the optimization method (e.g. running time etc)
    diagnostic.iterations = 0;
    diagnostic.residualHistory = [];
    diagnostic.runningTime = 0;
    
    % transfer the pts in all views together with the list of initial
    % homographies to
    % a common global coordinate system for numerical stability
    [sceneDataJointlyNormalised, ~, T, Tp] = ...
          compute_global_normalisation_transform(sceneData,listOfInitialH);
   
        
    % Compute DLT Estimates because we will need it as a seed for FNS
    % estimate
    [listOfNormalisedInitialH, ~] = ...
          compute_dlt_estimates_and_covariances_aml(...
                                        sceneDataJointlyNormalised, T, Tp);

    % Compute FNS Estimate with FNS Covariances 
    if(useFNSCovariances)
        [listOfNormalisedInitialH, ~] = ...
          compute_fns_estimates_and_covariances(...
                   sceneDataJointlyNormalised, 5,listOfNormalisedInitialH);
    end    
     
    % Add some outliers to the normalised scene data
    sceneDataJointlyNormalised = add_outliers_to_sceneData(...
        sceneDataJointlyNormalised ,0.3,[-1,1]);
    
    % check how many data points there are in the scene
    % and convert data points to homogenous coordinates
    totalNumberOfDataPoints = 0;
    for i = 1:length(sceneDataJointlyNormalised)
        numberOfPointsOnPlaneI = ...
                        length(sceneDataJointlyNormalised(i).ptsInViewOne);
        totalNumberOfDataPoints = ...
                          totalNumberOfDataPoints + numberOfPointsOnPlaneI;        
        sceneDataJointlyNormalised(i).ptsInViewOne = ...
               [sceneDataJointlyNormalised(i).ptsInViewOne ; ...
               ones(1,length(sceneDataJointlyNormalised(i).ptsInViewOne))];
        sceneDataJointlyNormalised(i).ptsInViewTwo = ...
               [sceneDataJointlyNormalised(i).ptsInViewTwo ; ...
               ones(1,length(sceneDataJointlyNormalised(i).ptsInViewTwo))];
    end
    
 

	% setup matlabs optimisation structs
    options = optimset('MaxIter',10000,'MaxFunEvals',100000, ...
        'FinDiffType ','central','TolFun',1e-10,'TolX',1e-10,...
        'Diagnostics','off','Display','iter','Algorithm',...
        'levenberg-marquardt','Jacobian','on','DerivativeCheck','off',...
        'OutputFcn', @myoutput);
    f = @(bestEta)cost_function_vector_sampsonaml(...
              bestEta,sceneDataJointlyNormalised,totalNumberOfDataPoints, ...
              threshold);
    
    if (useChojnackiInitialisation)
      % using globablly normalised data as input to initialisation
      eta = initialise_eta_chojnacki(listOfNormalisedInitialH);
    else
      % using globablly normalised data as input to initialisation
      eta = initialise_eta_peichen_true(listOfNormalisedInitialH);  
    end
    
    % keeps track of the residuals for each iterations of the optimisation
    residualHistory = [];
    
    fprintf('\n Running Robust Sampson AML Multiview Estimate \n')
    tic
	% Levenberg-Marquardt non-liner least squares optimisation
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

	% Convert the latent variables encoding all the 
    % homographies into individual homography matrices
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
%   consistent homographies into individual homography matrices.
%
%   Parameters:
%
%      eta          - vector of latent variables encoding the consistent
%                     homographies
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


%   Function: cost_function_vector_robustsampsonaml
%
%   This function computes the Sampson distance residual given a collection
%   of corresponding points associated with 'I' different homographies, 
%   and a vector of latent variables parametrising the 'I' different
%   homographies so that all inter-homography constraints are satisfied. 
%   
%
%   Parameters:
%
%      eta                           - vector of latent variables encoding 
%                                      the consistent homographies
%
%      sceneDataJointlyNormalised    - a structure containing corresponding
%                                      points for each homography all 
%                                      jointly normalised in a global 
%                                      coordinate system.
%
%     totalNumberOfDataPoints        - the sum of all corresponding points
%                                      for every homography
%										
%
%   Returns: (1) a vector containing the robust Sampson distance based 
%            residuals for all corresponding points  and (2) 
%            the jacobian matrix of the cost function
% 
function [residualVector, jacobianMatrix] = ...
                  cost_function_vector_sampsonaml(...
                    eta,sceneDataJointlyNormalised,totalNumberOfDataPoints,...
                    threshold)

residualVector = zeros(totalNumberOfDataPoints,1);
jacobianMatrix = zeros(totalNumberOfDataPoints,length(eta));

% there are 4*I+12 parameters for I planes
numberOfPlanes = (length(eta) - 12) / 4;
% extract first latent variable
a = eta(1:9);
% convert length-9 vector into 3x3 matrix
A = reshape(a,3,3);
% extact second latent variable
b = eta(10:12);

% note that A is not a true rotation matrix and b is not a true translation
% because the calibration matrix of the cameras is unknown. Refer to
% equations (1), (2) and (3) in our paper for more details. 

index = 13;
% Extract remaining latent variables
for i = 1:numberOfPlanes
    vi = eta(index:index+2);
    wi = eta(index+3);
    index = index + 4;
    if (i==1)
        CompositeHomographies = (wi * A + b * vi'); 
    else
        CompositeHomographies = ...
                            horzcat(CompositeHomographies,(wi*A + b* vi')); 
    end
    
    %alternative forumlation NB this saves us from having to do
    %reshape all the matrices.. but it is less transparent!
    %alt = vi*a' + kron(vn_i, b)
    listVn{i} = vi;
    listWi{i} = wi;    
end

% Vectorize the composite homography matrices
vectorOfCompositeHomographies = ...
                      reshape(CompositeHomographies,9 * numberOfPlanes,1);

dataPointCounter = 1;
homographyCounter = 1;
combinedJacobian = [];
% each homography matrix consists of 9 parameters
for k = 1:9:numberOfPlanes*9
    
    % choose the correct length 9 subset
    theta_i = vectorOfCompositeHomographies(k:k+8) ;
    
    % all data points associated with the k'th homography in view 1
    all_m = sceneDataJointlyNormalised(homographyCounter).ptsInViewOne;
    % all data points associated with the k'th homography in view 2
    all_mprime = sceneDataJointlyNormalised(homographyCounter).ptsInViewTwo;
    
    % loop over all corresponding points
    for iHomography = 1:length(all_m)
        % compute data carrier matrix (refer to section on "Background on
        % homography Model" for the origin of these equations
        m = all_m(:,iHomography); 
        mprime = all_mprime(:,iHomography);
        anti_mprime = [0 -mprime(3) mprime(2) ;
                       mprime(3) 0 -mprime(1) ;
                      -mprime(2) mprime(1) 0];
        Uz = kron(-m, anti_mprime);
        e1 = [1 ; 0 ; 0];
        e2 = [0 ; 1 ; 0];
        I = [e1 e2];
        Vz = Uz*I;
        
        % Compute Mahalanobis matrix
        anti_e1 = [0 -e1(3) e1(2) ; e1(3) 0 -e1(1) ; -e1(2) e1(1) 0];
        anti_e2 = [0 -e2(3) e2(2) ; e2(3) 0 -e2(1) ; -e2(2) e2(1) 0];
        derivVz = -[vec(kron(e1,anti_mprime)*I), ...
                    vec(kron(e2,anti_mprime)*I), ...
                    vec(kron(m,anti_e1)*I), ...
                    vec(kron(m,anti_e2)*I)];
        % this represents the covariance associated with the original
        % data points (i.e. the uncertainty in the x and y coordinates)
        % we assume that it is isotropic
        % due to the scale invariance of our AML objective function, the
        % scale of the covariance matrix is not important
        covariancekj = eye(4,4);
        Bkj = derivVz * covariancekj * derivVz';
        Sigmakj = kron(eye(2),theta_i') * Bkj * kron(eye(2),theta_i);
        
        
        % threshold for humber norm
        %threshold = 0.005;
        %threshold = Inf;
        sqrtSampson = sqrt(theta_i' * Vz * (Sigmakj \ Vz') * theta_i);
        residualVector(dataPointCounter)  = huberD(sqrtSampson,threshold);
%         if (abs(residual) < threshold)
%             
%         residualVector(dataPointCounter) = residual;
%         1;
%         else
%         residualVector(dataPointCounter) = ...
%                           sqrt((2*threshold)*abs(residual) - threshold^2);
%         1;
%         end
        
        %residualVector(dataPointCounter) = ...
        %                   sqrt(theta_i' * Vz * (Sigmakj \ Vz') * theta_i);
        
        
        % the following matrices will be used to compute the jacobian matrix
        % the names of the matrices and variables try to follow the
        % notation in the paper
        M =  Vz * (Sigmakj \ Vz');
        smallsigmakj = (Sigmakj \ Vz')*theta_i;
        N =  kron(smallsigmakj',eye(9)) * Bkj * kron(smallsigmakj,eye(9));
        X = M-N;
%         daf = listWi{homographyCounter} * ...
%                         residualVector(dataPointCounter)^-1 * theta_i' * X;
%         dbf = residualVector(dataPointCounter)^-1 * ...
%                      theta_i' * X * kron(listVn{homographyCounter},eye(3));
%         dvnf = residualVector(dataPointCounter)^-1 * ...
%                                             theta_i' * X *(kron(eye(3),b));
%         dvif = residualVector(dataPointCounter)^-1 * ...
%                                                          theta_i' * X * a ;
        daf = listWi{homographyCounter} * ...
                        sqrtSampson^-1 * theta_i' * X;
        dbf = sqrtSampson^-1 * ...
                     theta_i' * X * kron(listVn{homographyCounter},eye(3));
        dvnf = sqrtSampson^-1 * ...
                                            theta_i' * X *(kron(eye(3),b));
        dvif = sqrtSampson^-1 * ...
                                                         theta_i' * X * a ;
        
        % Combine all partial derivatives into one jacobian matrix
        myJacobian = horzcat(daf,dbf);
        for iHomography = 1:numberOfPlanes
            if (iHomography == homographyCounter)
                myJacobian = horzcat(myJacobian,dvnf);
                myJacobian = horzcat(myJacobian,dvif);
            else
                % the derivative is zero for latent variables that are not
                % associated with current homography
                myJacobian = horzcat(myJacobian,zeros(1,3));
                myJacobian = horzcat(myJacobian,zeros(1,1));
            end
        end
        
        % include the partial derivative w.r.t to huber norm composition 
        dPrime = huberDPrime( sqrtSampson,threshold);
        myJacobian = myJacobian*dPrime;
        length(combinedJacobian);
            length(myJacobian);
            combinedJacobian = vertcat(combinedJacobian,myJacobian);

                    
        dataPointCounter = dataPointCounter + 1;
    end       
    homographyCounter = homographyCounter + 1;
end

residualVector;

% return the jacobian
jacobianMatrix = combinedJacobian;

1;
end

function val = huberC(r,threshold)
if (abs(r) <= threshold)
    val = r^2;
else
    val = 2*threshold*abs(r) - threshold^2;
end
end

function val = huberCPrime(r, threshold)
    if (abs(r)) <= threshold
        val = 2*r;
    elseif (threshold < r)
        val = 2*threshold;
    elseif (r < -threshold)
        val = -2*threshold;
    end     
end

function val = huberD(r,threshold)
    %val = sqrt(huberC(r,threshold));
    if (abs(r)) <= threshold
        val = r;
    elseif (r > threshold )        
        val = sqrt(2*threshold*r-threshold^2);
    elseif (r < -threshold)
        val = -sqrt(-2*threshold*r-threshold^2);
    end 
end

function val = huberDPrime(r,threshold)
    if (abs(r)) <= threshold
        val = 1;
    elseif (r > threshold)
        %val = threshold*(2*threshold*r-threshold^2)^-0.5;
        val = threshold / sqrt(2*threshold*r-threshold^2);
    elseif (r < -threshold)
        val = threshold*(-2*threshold*r-threshold^2)^-0.5;
    end 
end


