%   Function: compute_aml_estimates
%
%   Given a collection of homographies estimated using DLT or FNS together 
%   with covariance matrices, this function 'upgrades' the homographies so
%   that they satisfy all possible two-view consistency   constraints that 
%   arise due to the rigidity of the scene. For example, two homographies
%   corresponding to two facades of a building can be made fully consistent
%   because the building is rigid. However,  if you somehow model the side
%   views of two cars with two different homographies, and these
%   cars change their relative positions then these homographies cannot be
%   made fully consistent because the assumption of rigidity is violated.
%
%   This algorithm originally appeared in
%
%   Chojnacki, W.; Szpak, Z.L.; Brooks, M.J.; van den Hengel, A.; ,
%   "Multiple Homography Estimation with Full Consistency Constraints,"
%   Digital Image Computing: Techniques and Applications (DICTA), 2010
%   International Conference on , vol., no., pp.480-485, 1-3 Dec. 2010
%   doi: 10.1109/DICTA.2010.87
%
%   A further evaluation of the method appeared in
%
%   Z. L. Szpak, W. Chojnacki, A. Eriksson, and A. van den Hengel. 
%   Sampson distance based joint estimation of multiple homographies with
%   uncalibrated cameras. 
%   Comput. Vis. Image Underst., 125:200-213, 2014. 
%   http://dx.doi.org/10.1016/j.cviu.2014.04.008
%
%
%   At the time of publication this method is the fastest method for
%   enforcing FULL homography consistency constraints between two views.
%   For small noise levels, the performance of this method is almost 
%   indistinguishable from a gold standard variant of bundle adjustment 
%   that enforces full consistency constraints.
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
%      useChojnackiInitialisation  - simple flag variable, toggling between
%                                    using a Chojnacki initialisation to
%                                    initialise the required latent 
%                                    variables or Chen's initialisation
%                                    method. Chojnacki's method works
%									 when there are two or more 
%                                     homographies in the scene, while
%									 Chen's requires at least three 
%                                    homographies.
%
%      useFNSCovariances           - simple flag variable, toggling between
%                                    using FNS covariances or DLT 
%                                    covariances.   
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
                      compute_aml_estimates(sceneData, ...
                                        useChojnackiInitialisation,...
                                         useFNSCovariances)
    
    % data structure used to keep track of diagnostic information relating
    % to the performance of the optimization method (e.g. running time etc)
    diagnostic.iterations = 0;
    diagnostic.residualHistory = [];
    diagnostic.runningTime = 0;
    
    
    % transfer the pts in all views and the list of initial homographies to
    % a common global coordinate system for numerical stability
    [sceneDataJointlyNormalised, ~, T, Tp] = ...
         compute_global_normalisation_transform(sceneData);

    
    % Compute DLT Estimates because we will need it as a seed for FNS
    % estimate
    [listOfNormalisedInitialH, listOfCovariancesOfH] = ...
                     compute_dlt_estimates_and_covariances_aml(...
                                        sceneDataJointlyNormalised);

    % Compute FNS Estimate with FNS Covariances and use those for AML 
    if(useFNSCovariances)
        maxFnsIterations = 5;
        [listOfNormalisedInitialH, listOfCovariancesOfH] = ...
            compute_fns_estimates_and_covariances(...
                   sceneDataJointlyNormalised, maxFnsIterations, ...
                                                 listOfNormalisedInitialH);
    end    
     
    
    X = cell2mat(listOfNormalisedInitialH);
	
    % long vector of all initial homography matrices
    x_all = reshape(X,length(listOfNormalisedInitialH)*9,1);
    
	% setup matlabs optimisation structs
    options = optimset('FinDiffType ','central','Diagnostics','on',...
                    'Display','iter','Algorithm','levenberg-marquardt',...
                                            'TolFun',1e-10,'TolX',1e-18,...
                                'Jacobian','on','DerivativeCheck','off',...
                                'OutputFcn', @myoutput,...
                                'MaxIter',10000);
    f = @(bestEta)cost_function_vector_aml(...
                                       bestEta,x_all,listOfCovariancesOfH);
    
    if (useChojnackiInitialisation)
      % using globablly normalised data as input to initialisation
      eta = initialise_eta_chojnacki(listOfNormalisedInitialH);
    else
      % using globablly normalised data as input to initialisation
      eta = initialise_eta_peichen_true(listOfNormalisedInitialH);  
    end
    
    % keeps track of the residuals for each iterations of the optimisation
    residualHistory = [];
    
    fprintf('\n Running AML Multiview Estimate \n')
    tic % start timer
	% Levenberg-Marquardt non-linear least squares optimisation
    [etaFinal,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                            lsqnonlin(f,eta,[],[],options);
    runningTime = toc; % end timer                                  
    % Nested function whose sole purpose is to keep track of the
    % AML cost function at each iteration. The values will be used later
    % to analyse the convergence properties of our algorithm
    function stop = myoutput(~,optimvalues,state)
        stop = false;
        if state == 'iter'
          residualHistory = [residualHistory; optimvalues.resnorm];
        end
    end
                                    
    
	% Convert the latent variables encoding all the homographies
    % into individual homography matrices
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


%   Function: cost_function_vector_aml
%
%   This function computes the Sampson distance residual given a collection
%   of corresponding points associated with 'I' different homographies, 
%   and a vector of latent variables parametrising the 'I' different
%   homographies so that all inter-homography constraints are satisfied. 
%
%   The structure of the eta parameter (the latent variables) is slightly
%   different from the notation used in the paper. 
%   
%   The eta vector is : [a, b, v_1,w_1, ... , v_n, w_n], whereas in the
%   paper we structure it as [w, b, v_1, ... , v_n, w_1, ...,w_n]; 
%   Keep in mind that this has implications for the overall structure
%   of the Jacobian matrix.
%
%
%   Parameters:
%
%      eta 	   			     - vector of latent variables encoding the 
%                              consistent homographies
%
%      x_all  			     - a vector consisting of the original 
%                              'inconsistent' homographies
%                              all concatenated together
%
%      listOfCovariancesOfH  - cell array containing homography covariances 
%										
%
%   Returns: A matrix containing the cost function values for each
%            homography and a jacobian matrix
% 
function [matrixOfCostOfEachPlane,jacobianMatrix] = ...
                   cost_function_vector_aml(eta,x_all,listOfCovariancesOfH)
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
% Extract the plane normals and distance from origin
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

% Vectorize the homography matrices
vectorOfCompositeHomographies = ...
                    reshape(CompositeHomographies,9 * numberOfPlanes,1);
count = 1;
% each homography matrix consists of 9 parameters
for k = 1:9:numberOfPlanes*9
    
    % choose the correct length 9 subset
    theta_i = vectorOfCompositeHomographies(k:k+8) ;
    
    % alternative formula....
    %altTheta_i =  distancesVi{count} * a + kron(normalsVn{count},b);
    
    x_i = x_all(k:k+8);
    
    covariance_i = listOfCovariancesOfH{count};    

    % Symmetric projection matrix
    Pox_i = eye(9,9) - (x_i * x_i')/ (norm(x_i,2)^2) ;    

    Lambda = Pox_i * covariance_i * Pox_i;    
    
	% rank-8 constrained pseudo-inverse
    [U, S, V] = svd(Lambda);
    for n = 1:8
        S(n,n) = 1/(S(n,n))  ;        
    end
    S(9,9) = 0;
    
    Ai =  V * S *  U';    

    [U, S, V] = svd(Ai);
    for n = 1:9
        S(n,n) = sqrt(S(n,n))  ;
    end
    Bi = U*S*V';
    
    %disp('Cost of Fi')
    fi = (norm(theta_i,2))^-1 * (Bi*theta_i);
    ans = fi' * fi;
    %ans2 = norm(fi)
    
    % The sum of squares should not be formed explicitly,
    % instead we will need to return a vector of function
    % values. This is what matlabs optimisation requires.
    if (count == 1)
        %disp('Matrix of Costs for fi one')
        matrixOfCosts =  fi;
    else
        %disp('Matrix of Costs for fi two')
        matrixOfCosts = horzcat(matrixOfCosts,fi);
    end    
    %disp('This is the matrix of costs')
    %matrixOfCosts
    
    %_-_-_-_-_-_-_-_-_-building the jacobian matrix _-_-_-_-_-_-_-_-_-_-_- 
    Potheta_i = eye(9,9) - (theta_i * theta_i')/ (norm(theta_i,2)^2) ;
    a;
    b;
    Potheta_i;
    listVn{count};
    Bi;    
    % partial derivatives
    daf = listWi{count} * (1/norm(theta_i,2)) * Bi * Potheta_i  ;
    dbf = ...
       (1/norm(theta_i,2)) * Bi * Potheta_i * (kron(listVn{count},eye(3)));
    dvnf = (1/norm(theta_i,2)) * Bi * Potheta_i*(kron(eye(3),b));
    dvif = (1/norm(theta_i,2)) * Bi * Potheta_i * a   ;
    
    myJacobian = horzcat(daf,dbf);
    for j = 1:numberOfPlanes
        if (j == count)
            myJacobian = horzcat(myJacobian,dvnf);
            myJacobian = horzcat(myJacobian,dvif);
        else
            myJacobian = horzcat(myJacobian,zeros(9,3));
            myJacobian = horzcat(myJacobian,zeros(9,1));
        end
    end
    
    if (count == 1)
        combinedJacobian = myJacobian;
    else
        length(combinedJacobian);
        length(myJacobian);
        combinedJacobian = vertcat(combinedJacobian,myJacobian);
    end    
    %_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ 
	
    count = count+1;
end

% return the vector of cost for each plane
matrixOfCostOfEachPlane = matrixOfCosts;
% return the jacobian
jacobianMatrix = combinedJacobian;
end

