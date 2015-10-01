%   Function: compute_sampsonaml_estimates
%
%   Given a collection of homographies estimated using DLT or FNS together with covariance matrices,
%   this function 'upgrades' the homographies so that they satisfy all possible two-view consistency
%   constraints that arise due to the rigidity of the scene. For example, two homographies corresponding
%   to two facades of a building can be made fully consistent because the building is rigid. However,
%   if you somehow model the side views of two cars with two different homographies, and these
%   cars change their relative positions then these homographies cannot be made fully consistent. 
%
%   Parameters:
%
%      sceneData 	   			   - a struct containing noisy corresponding points between two views
%
%      listOfInitialH  			   - initial homography estimates, either DLT or FNS based
%
%      useChojnackiInitialisation  - simple flag variable, toggling between using a Chojnacki
%									 initialisation to initialise the required latent variables
%								     or Chen's initialisation method. Chojnacki's method works
%									 when there are two or more homographies in the scene, while
%									 Chen's requires at least three homographies.
%
%      useFNSCovariances           - simple flag variable, toggling between using FNS covariances
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

function listOfEstimatedH = compute_sampsonaml_estimates(sceneData, listOfInitialH, useChojnackiInitialisation, useFNSCovariances )

    
    % transfer the pts in all views and the list of initial homographies to
    % a common global coordinate system for numerical stability
    [sceneDataJointlyNormalised ignore_me T Tp] = compute_global_normalisation_transform(sceneData,listOfInitialH);
   
    listOfCovariancesOfH = cell(0,length(listOfInitialH));
    
    % Compute DLT Estimates because we will need it as a seed for FNS
    % estimate
    [listOfNormalisedInitialH listOfCovariancesOfH] = compute_dlt_estimates_and_covariances_aml(sceneDataJointlyNormalised, T, Tp);

    % Compute FNS Estimate with FNS Covariances and use those for AML 
    if(useFNSCovariances)
        [listOfNormalisedInitialH listOfCovariancesOfH] = compute_fns_estimates_and_covariances(sceneDataJointlyNormalised, 5,listOfNormalisedInitialH);
    end    
     
    totalNumberOfDataPoints = 0;
    % check how many data points there are in the scene
    % and convert data points to homogenous coordinates
    for i = 1:length(sceneDataJointlyNormalised)
        numberOfPointsOnPlaneI = length(sceneDataJointlyNormalised(i).ptsInViewOne);
        totalNumberOfDataPoints = totalNumberOfDataPoints + numberOfPointsOnPlaneI;        
        sceneDataJointlyNormalised(i).ptsInViewOne = [sceneDataJointlyNormalised(i).ptsInViewOne ; ones(1,length(sceneDataJointlyNormalised(i).ptsInViewOne))];
        sceneDataJointlyNormalised(i).ptsInViewTwo = [sceneDataJointlyNormalised(i).ptsInViewTwo ; ones(1,length(sceneDataJointlyNormalised(i).ptsInViewTwo))];
    end
    
 

	% setup matlabs optimisation structs
    options = optimset('MaxIter',10000,'MaxFunEvals',100000, 'FinDiffType ','central','TolFun',1e-10,'TolX',1e-10,'Diagnostics','off','Display','iter','Algorithm','levenberg-marquardt','Jacobian','on','DerivativeCheck','off');
    f = @(bestEta)cost_function_vector_sampsonaml(bestEta,sceneDataJointlyNormalised,totalNumberOfDataPoints);
    
    if (useChojnackiInitialisation)
      % using globablly normalised data as input to initialisation
      eta = initialise_eta_chojnacki(listOfNormalisedInitialH);
    else
      % using globablly normalised data as input to initialisation
      eta = initialise_eta_peichen_true(listOfNormalisedInitialH);  
    end
    
    %debugging derivative check
   % eta = ones(length(eta),1);
%     eta = [
%           1.02707642036823;
%           0.93003492368847;
%           1.02400881344689;
%          0.997901890210339;
%           1.02308305572916;
%           1.01992157865388;
%            0.9705316415329;
%          0.949670353919302;
%            1.0537704779361;
%          0.995647563905817;
%          0.903196177839156;
%           1.09755301964135;
%           1.00012453864907;
%           1.01109343907794;
%          0.989089865597024;
%          0.999633737784883;
%          0.983976844601953;
%           1.02922388400107;
%          0.987203292801941;
%          0.999281677090805];
    
    fprintf('\n Running AML Multiview Estimate \n')
	% Levenberg-Marquardt non-liner least squares optimisation
    [etaFinal,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(f,eta,[],[],options);
    
    
    % [residualVector jacobianMatrix]= cost_function_vector_sampsonaml(etaFinal,sceneDataJointlyNormalised,totalNumberOfDataPoints)
     %J = finjac(f,residualVector,eta,1e-10)
     
   %  J = finjac(f,residual,etaFinal,1e-10);
    
	% Convert the latent variables encoding all the homographies into individual homography matrices
    listOfEstimatedH = convertEtaToHomographyCellArray(etaFinal);
    
    
    for i = 1:length(listOfEstimatedH)
        listOfEstimatedH{i} =  Tp \ listOfEstimatedH{i} * T;
    end
    
end


%   Function: convertEtaToHomographyCellArray
%
%   Converts a vector of latent variables encoding a collection of fully consistent homographies
%   into individual homography matrices.
%
%   Parameters:
%
%      eta 	   			     - vector of latent variables encoding the consistent homographies
%
%										
%
%   Returns: A cell array containing fully consistent homography matrices
% 
function listOfOptimisedValidHomographies = convertEtaToHomographyCellArray(eta)
% HARD CODED FOR 4 PLANES AT THE MOMENT FIX ME FIX ME!
numberOfPlanes = (length(eta) - 12) / 4;
listOfOptimisedValidHomographies = cell(1,numberOfPlanes);

a = eta(1:9);
% vertical concatenating the rows to form 3 into 3
% rotation matrix
A = reshape(a,3,3);
% extract the translation of the camera
b = eta(10:12);

index = 13;
% Extract the plane normals and distance from origin
for i = 1:numberOfPlanes
    vn_i = eta(index:index+2);
    vi = eta(index+3);
    index = index + 4;
    listOfOptimisedValidHomographies{i} = (vi * A + b * vn_i'); 
    
end

end


%   Function: cost_function_vector_aml
%
%   Given a collection of homographies estimated using DLT or FNS together with covariance matrices,
%
%   Parameters:
%
%      eta 	   			     - vector of latent variables encoding the consistent homographies
%
%      x_all  			     - a vector consisting of the original 'inconsistent' homographies
%                              all concatenated together
%
%      listOfCovariancesOfH  - cell array containing homography covariances 
%										
%
%   Returns: A matrix containing the cost function values for each homography and a jacobian matrix
% 
function [residualVector, jacobianMatrix] = cost_function_vector_sampsonaml(eta,sceneDataJointlyNormalised,totalNumberOfDataPoints)

residualVector = zeros(totalNumberOfDataPoints,1);
jacobianMatrix = zeros(totalNumberOfDataPoints,length(eta));

numberOfPlanes = (length(eta) - 12) / 4;
a = eta(1:9);
% vertical concatenating the rows to form 3 into 3
% 'rotation matrix'
A = reshape(a,3,3);
% extract the 'translation' of the camera
b = eta(10:12);

% note that A is not a true rotation matrix and b is not a true translation
% because the calibration matrix of the cameras is unknown. So A is actually the
% product of a rotation matrix and a calibrtation matrix and similiarily b is also
% a product of the translation and calibration matrix.

% from this point onward, the latent variables encode the planes
index = 13;
% Extract the plane normals and distance from origin
for i = 1:numberOfPlanes
    vn_i = eta(index:index+2);
    vi = eta(index+3);
    index = index + 4;
    if (i==1)
        THETA = (vi * A + b * vn_i'); 
    else
        THETA = horzcat(THETA,(vi*A + b* vn_i')); 
    end
    
    %alternative forumlation NB this saves us from having to do
    %reshape all the matrices.. but it is less transparent!
    %alt = vi*a' + kron(vn_i, b)
    normalsVn{i} = vn_i;
    distancesVi{i} = vi;    
end

% Vectorize the homography matrices
theta_all = reshape(THETA,9 * numberOfPlanes,1);
count = 1;
dataPointCounter = 1;
homographyCounter = 1;
jacobianCounter = 1;
combinedJacobian = [];
for k = 1:9:numberOfPlanes*9
    
    % choose the correct length 9 subset
    theta_i = theta_all(k:k+8) ;
    
    all_m = sceneDataJointlyNormalised(homographyCounter).ptsInViewOne;
    all_mprime = sceneDataJointlyNormalised(homographyCounter).ptsInViewTwo;
    M = zeros(9,9);
    N = zeros(9,9);
    X = zeros(9,9);
    for j = 1:length(all_m)
        % compute data carrier matrix
        m = all_m(:,j); 
        mprime = all_mprime(:,j);
        anti_mprime = [0 -mprime(3) mprime(2) ; mprime(3) 0 -mprime(1) ; -mprime(2) mprime(1) 0];
        Uz = kron(-m, anti_mprime);
        e1 = [1 ; 0 ; 0];
        e2 = [0 ; 1 ; 0];
        I = [e1 e2];
        Vz = Uz*I;
        
        % compute mahalanobis matrix
        anti_e1 = [0 -e1(3) e1(2) ; e1(3) 0 -e1(1) ; -e1(2) e1(1) 0];
        anti_e2 = [0 -e2(3) e2(2) ; e2(3) 0 -e2(1) ; -e2(2) e2(1) 0];
        derivVz = -[vec(kron(e1,anti_mprime)*I), vec(kron(e2,anti_mprime)*I) , vec(kron(m,anti_e1)*I), vec(kron(m,anti_e2)*I)];
        covariancekj = eye(4,4);
        Bkj = derivVz * covariancekj * derivVz';
        Sigmakj = kron(eye(2),theta_i') * Bkj * kron(eye(2),theta_i);
     

        1;
        
        %res = theta_i' * Vz * inv(Sigmakj) * Vz' * theta_i;
        % numerically more stable
        %res = theta_i' * Vz * (Sigmakj \ Vz') * theta_i;
        1;
        residualVector(dataPointCounter) = sqrt(theta_i' * Vz * (Sigmakj \ Vz') * theta_i);
        
        % these matrices will be used to compute the jacobian matrix
        M =  Vz * (Sigmakj \ Vz');
        smallsigmakj = (Sigmakj \ Vz')*theta_i;

%         M =  Vz * inv(Sigmakj) * Vz';
%         smallsigmakj = inv(Sigmakj)* Vz'*theta_i;
        N =  kron(smallsigmakj',eye(9)) * Bkj * kron(smallsigmakj,eye(9));
        X = M-N;
         daf = distancesVi{homographyCounter} * residualVector(dataPointCounter)^-1 * theta_i' * X;
        dbf = residualVector(dataPointCounter)^-1 * theta_i' * X * kron(normalsVn{homographyCounter},eye(3));
        dvnf = residualVector(dataPointCounter)^-1 * theta_i' * X *(kron(eye(3),b));
        dvif = residualVector(dataPointCounter)^-1 * theta_i' * X * a ;


        
        % combine all partial derivatives
        myJacobian = horzcat(daf,dbf);
        for j = 1:numberOfPlanes
            if (j == homographyCounter)
                myJacobian = horzcat(myJacobian,dvnf);
                myJacobian = horzcat(myJacobian,dvif);
            else
                myJacobian = horzcat(myJacobian,zeros(1,3));
                myJacobian = horzcat(myJacobian,zeros(1,1));
            end
        end
        
            length(combinedJacobian);
            length(myJacobian);
            combinedJacobian = vertcat(combinedJacobian,myJacobian);

        
        %jacobianCounter = jacobianCounter + 1;
                
        dataPointCounter = dataPointCounter + 1;
    end
    
%     X =  M - N;
%     % compute jacobian matrix
%     for dummyCounter = 1:length(all_m)
%         daf = distancesVi{homographyCounter} * residualVector(jacobianCounter)^-1 * theta_i' * X;
%         dbf = residualVector(jacobianCounter)^-1 * theta_i' * X * kron(normalsVn{homographyCounter},eye(3));
%         dvnf = residualVector(jacobianCounter)^-1 * theta_i' * X *(kron(eye(3),b));
%         dvif = residualVector(jacobianCounter)^-1 * theta_i' * X * a ;
% 
% 
%         
%         % combine all partial derivatives
%         myJacobian = horzcat(daf,dbf);
%         for j = 1:numberOfPlanes
%             if (j == homographyCounter)
%                 myJacobian = horzcat(myJacobian,dvnf);
%                 myJacobian = horzcat(myJacobian,dvif);
%             else
%                 myJacobian = horzcat(myJacobian,zeros(1,3));
%                 myJacobian = horzcat(myJacobian,zeros(1,1));
%             end
%         end
%         
% %         if (homographyCounter == 1)
% %             combinedJacobian = myJacobian;
% %         else
%             length(combinedJacobian);
%             length(myJacobian);
%             combinedJacobian = vertcat(combinedJacobian,myJacobian);
%         %end
%         
%         jacobianCounter = jacobianCounter + 1;
%         
%     end
%     

    
    homographyCounter = homographyCounter + 1;
    
    %     % partial derivatives
%     %daf = distancesVi{count} * (norm(theta_i,2)^-1) * Bi * Potheta_i  ;
%     daf = distancesVi{count} * (1/norm(theta_i,2)) * Bi * Potheta_i  ;
%     %dbf = (norm(theta_i,2)^-1) * Bi * Potheta_i * (kron(normalsVn{count},eye(3)));
%     dbf = (1/norm(theta_i,2)) * Bi * Potheta_i * (kron(normalsVn{count},eye(3)));
%     %dvnf = (norm(theta_i,2)^-1) * Bi * Potheta_i*(kron(eye(3),b));
%     dvnf = (1/norm(theta_i,2)) * Bi * Potheta_i*(kron(eye(3),b));
%     %dvif = (norm(theta_i,2)^-1) * Bi * Potheta_i * a   ;
%     dvif = (1/norm(theta_i,2)) * Bi * Potheta_i * a   ;
%     
%     myJacobian = horzcat(daf,dbf);
%     for j = 1:numberOfPlanes
%         if (j == count)
%             myJacobian = horzcat(myJacobian,dvnf);
%             myJacobian = horzcat(myJacobian,dvif);
%         else
%             myJacobian = horzcat(myJacobian,zeros(9,3));
%             myJacobian = horzcat(myJacobian,zeros(9,1));
%         end
%     end
%     
    
    1;
%     
%     
%     % alternative formula....
%     %altTheta_i =  distancesVi{count} * a + kron(normalsVn{count},b);
%     
%     x_i = x_all(k:k+8);
%     
%     covariance_i = listOfCovariancesOfH{count};    
% 
%     % symmetric projection matrix
%     Pox_i = eye(9,9) - (x_i * x_i')/ (norm(x_i,2)^2) ;    
% 
%     Lambda = Pox_i * covariance_i * Pox_i;    
%     
% 	% rank-8 constrained pseudo-inverse
%     [U S V] = svd(Lambda);
%       for n = 1:8
%         S(n,n) = 1/(S(n,n))  ;        
%     end
%     S(9,9) = 0;
%     
%     Ai =  V * S *  U';    
% 
%     [U S V] = svd(Ai);
%     for n = 1:9
%         S(n,n) = sqrt(S(n,n))  ;
%     end
%     Bi = U*S*V';
%     
%     %disp('Cost of Fi')
%     fi = (norm(theta_i,2))^-1 * (Bi*theta_i);
%     ans = fi' * fi;
%     %ans2 = norm(fi)
%     
%     % The sum of squares should not be formed explicitly,
%     % instead we will need to return a vector of function
%     % values. This is what matlabs optimisation requires.
%     if (count == 1)
%         %disp('Matrix of Costs for fi one')
%         matrixOfCosts =  fi;
%     else
%         %disp('Matrix of Costs for fi two')
%         matrixOfCosts = horzcat(matrixOfCosts,fi);
%     end    
%     %disp('This is the matrix of costs')
%     %matrixOfCosts
%     
%     %////////////////////// building the jacobian matrix //////////////////////////
%     Potheta_i = eye(9,9) - (theta_i * theta_i')/ (norm(theta_i,2)^2) ;
%     %disp('a');
%     a;
%     %disp('b');
%     b;
%     % disp('Potheta_i');
%     Potheta_i;
%     %disp('normalsVn');
%     normalsVn{count};
%     %disp('Bi');
%     Bi;    
%     % partial derivatives
%     %daf = distancesVi{count} * (norm(theta_i,2)^-1) * Bi * Potheta_i  ;
%     daf = distancesVi{count} * (1/norm(theta_i,2)) * Bi * Potheta_i  ;
%     %dbf = (norm(theta_i,2)^-1) * Bi * Potheta_i * (kron(normalsVn{count},eye(3)));
%     dbf = (1/norm(theta_i,2)) * Bi * Potheta_i * (kron(normalsVn{count},eye(3)));
%     %dvnf = (norm(theta_i,2)^-1) * Bi * Potheta_i*(kron(eye(3),b));
%     dvnf = (1/norm(theta_i,2)) * Bi * Potheta_i*(kron(eye(3),b));
%     %dvif = (norm(theta_i,2)^-1) * Bi * Potheta_i * a   ;
%     dvif = (1/norm(theta_i,2)) * Bi * Potheta_i * a   ;
%     
%     myJacobian = horzcat(daf,dbf);
%     for j = 1:numberOfPlanes
%         if (j == count)
%             myJacobian = horzcat(myJacobian,dvnf);
%             myJacobian = horzcat(myJacobian,dvif);
%         else
%             myJacobian = horzcat(myJacobian,zeros(9,3));
%             myJacobian = horzcat(myJacobian,zeros(9,1));
%         end
%     end
%     
%     if (count == 1)
%         combinedJacobian = myJacobian;
%     else
%         length(combinedJacobian);
%         length(myJacobian);
%         combinedJacobian = vertcat(combinedJacobian,myJacobian);
%     end    
%     %/////////////////////////////////////////////////////////////////////////////
% 	
%     count = count+1;
end

residualVector;

% return the vector of cost for each plane
%matrixOfCostOfEachPlane = matrixOfCosts;
% return the jacobian
jacobianMatrix = combinedJacobian;


%     g = @(bestEta)cost_function_vector_sampsonaml_CLONE(bestEta,sceneDataJointlyNormalised,totalNumberOfDataPoints);
% J = finjac(g,residualVector,eta,1e-15)
% jacobianMatrix = J;
1;
end

function J = finjac(FUN,r,x,epsx)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lx=length(x);
J=zeros(length(r),lx);
epsx=epsx*ones(lx,1);
for k=1:lx
    dx=.25*epsx(k);
    xd=x;
    xd(k)=xd(k)+dx;
    rd=feval(FUN,xd);
%   ~~~~~~~~~~~~~~~~    
    J(:,k)=((rd-r)/dx);
end
end

% 
% function [residualVector jacobianMatrix] = cost_function_vector_sampsonaml_CLONE(eta,sceneDataJointlyNormalised,totalNumberOfDataPoints)
% 
% residualVector = zeros(totalNumberOfDataPoints,1);
% jacobianMatrix = zeros(totalNumberOfDataPoints,length(eta));
% 
% numberOfPlanes = (length(eta) - 12) / 4;
% a = eta(1:9);
% % vertical concatenating the rows to form 3 into 3
% % 'rotation matrix'
% A = reshape(a,3,3);
% % extract the 'translation' of the camera
% b = eta(10:12);
% 
% % note that A is not a true rotation matrix and b is not a true translation
% % because the calibration matrix of the cameras is unknown. So A is actually the
% % product of a rotation matrix and a calibrtation matrix and similiarily b is also
% % a product of the translation and calibration matrix.
% 
% % from this point onward, the latent variables encode the planes
% index = 13;
% % Extract the plane normals and distance from origin
% for i = 1:numberOfPlanes
%     vn_i = eta(index:index+2);
%     vi = eta(index+3);
%     index = index + 4;
%     if (i==1)
%         THETA = (vi * A + b * vn_i'); 
%     else
%         THETA = horzcat(THETA,(vi*A + b* vn_i')); 
%     end
%     
%     %alternative forumlation NB this saves us from having to do
%     %reshape all the matrices.. but it is less transparent!
%     %alt = vi*a' + kron(vn_i, b)
%     normalsVn{i} = vn_i;
%     distancesVi{i} = vi;    
% end
% 
% % Vectorize the homography matrices
% theta_all = reshape(THETA,9 * numberOfPlanes,1);
% count = 1;
% dataPointCounter = 1;
% homographyCounter = 1;
% jacobianCounter = 1;
% combinedJacobian = [];
% for k = 1:9:numberOfPlanes*9
%     
%     % choose the correct length 9 subset
%     theta_i = theta_all(k:k+8) ;
%     
%     all_m = sceneDataJointlyNormalised(homographyCounter).ptsInViewOne;
%     all_mprime = sceneDataJointlyNormalised(homographyCounter).ptsInViewTwo;
%     M = zeros(9,9);
%     N = zeros(9,9);
%     X = zeros(9,9);
%     for j = 1:length(all_m)
%         % compute data carrier matrix
%         m = all_m(:,j); 
%         mprime = all_mprime(:,j);
%         anti_mprime = [0 -mprime(3) mprime(2) ; mprime(3) 0 -mprime(1) ; -mprime(2) mprime(1) 0];
%         Uz = kron(-m, anti_mprime);
%         e1 = [1 ; 0 ; 0];
%         e2 = [0 ; 1 ; 0];
%         I = [e1 e2];
%         Vz = Uz*I;
%         
%         % compute mahalanobis matrix
%         anti_e1 = [0 -e1(3) e1(2) ; e1(3) 0 -e1(1) ; -e1(2) e1(1) 0];
%         anti_e2 = [0 -e2(3) e2(2) ; e2(3) 0 -e2(1) ; -e2(2) e2(1) 0];
%         derivVz = -[vec(kron(e1,anti_mprime)*I), vec(kron(e2,anti_mprime)*I) , vec(kron(m,anti_e1)*I), vec(kron(m,anti_e2)*I)];
%         covariancekj = eye(4,4);
%         Bkj = derivVz * covariancekj * derivVz';
%         Sigmakj = kron(eye(2),theta_i') * Bkj * kron(eye(2),theta_i);
%      
%         
%         %res = theta_i' * Vz * inv(Sigmakj) * Vz' * theta_i;
%         % numerically more stable
%         %res = theta_i' * Vz * (Sigmakj \ Vz') * theta_i;
%         1;
%         residualVector(dataPointCounter) = sqrt(theta_i' * Vz * (Sigmakj \ Vz') * theta_i);
%         
%         % these matrices will be used to compute the jacobian matrix
% %         M = M + Vz * (Sigmakj \ Vz');
% %         smallsigmakj = (Sigmakj \ Vz')*theta_i;
%         M = M + Vz * inv(Sigmakj) * Vz';
%         smallsigmakj = inv(Sigmakj)* Vz'*theta_i;
%         N = N + kron(smallsigmakj',eye(9)) * Bkj * kron(smallsigmakj,eye(9));
%                 
%         dataPointCounter = dataPointCounter + 1;
%     end
%     
%     X = M - N;
%     % compute jacobian matrix
%     for dummyCounter = 1:length(all_m)
%         daf = distancesVi{homographyCounter} * residualVector(jacobianCounter)^-1 * theta_i' * X;
%         dbf = residualVector(jacobianCounter)^-1 * theta_i' * X * kron(normalsVn{homographyCounter},eye(3));
%         dvnf = residualVector(jacobianCounter)^-1 * theta_i' * X *(kron(eye(3),b));
%         dvif = residualVector(jacobianCounter)^-1 * theta_i' * X * a ;
% 
% 
%         
%         % combine all partial derivatives
%         myJacobian = horzcat(daf,dbf);
%         for j = 1:numberOfPlanes
%             if (j == homographyCounter)
%                 myJacobian = horzcat(myJacobian,dvnf);
%                 myJacobian = horzcat(myJacobian,dvif);
%             else
%                 myJacobian = horzcat(myJacobian,zeros(1,3));
%                 myJacobian = horzcat(myJacobian,zeros(1,1));
%             end
%         end
%         
% %         if (homographyCounter == 1)
% %             combinedJacobian = myJacobian;
% %         else
%             length(combinedJacobian);
%             length(myJacobian);
%             combinedJacobian = vertcat(combinedJacobian,myJacobian);
%         %end
%         
%         jacobianCounter = jacobianCounter + 1;
%         
%     end
%     
% 
%     
%     homographyCounter = homographyCounter + 1;
%     
% 
% 
% residualVector;
% 
% % return the vector of cost for each plane
% %matrixOfCostOfEachPlane = matrixOfCosts;
% % return the jacobian
% jacobianMatrix = combinedJacobian;
% 
% end
% end