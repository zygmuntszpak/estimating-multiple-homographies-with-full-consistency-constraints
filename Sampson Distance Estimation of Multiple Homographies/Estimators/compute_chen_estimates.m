
% Computes an estimate of several homographies
% corresponding to several planes in a scene using the Weighted Alternate
% Least Squares method of Chen and Suter

% This code needs to be cleaned up as in its current form it includes
% several snippets of code that were attempts to improve the original
% WALS method

% 
function [listOfEstimatedH, diagnostic] = ...
     compute_chen_estimates(sceneData, listOfInitialH, useFNSCovariances )



    % transfer the pts in all views and the list of initial homographies to
    % a common coordinate system for numerical stability
    [sceneDataJointlyNormalised, ~, T, Tp] = ...
        compute_global_normalisation_transform(sceneData,listOfInitialH);
   
    listOfCovariancesOfH = cell(0,length(listOfInitialH));
    
    % Compute DLT Estimates because we will need it as a seed for FNS
    % estimate
    [listOfNormalisedInitialH, listOfCovariancesOfH] = ...
                compute_dlt_estimates_and_covariances_aml(...
                                        sceneDataJointlyNormalised);
    1;
    % Compute FNS Estimate with FNS Covariances and use those for AML 
    if(useFNSCovariances)
        [listOfNormalisedInitialH, listOfCovariancesOfH] = ...
                compute_fns_estimates_and_covariances(...
                  sceneDataJointlyNormalised, 5,listOfNormalisedInitialH);
    end    

 

    
    X = cell2mat(listOfNormalisedInitialH);
    % long vector of all initial homography matrices
    X_All = reshape(X,9,length(listOfNormalisedInitialH));
    
    
     [eta, listOfValidInitialH, A, b, vs, ws] = ...
                   initialise_eta_peichen_true(listOfNormalisedInitialH);
    
    % call optimization with eta    
    1;
    
    fprintf('\n Running Pei-Chen Estimate \n')
 
    %listOfEstimatedH = optimise_eta_peichen_algorithm_projective(A,b,vs,ws,X_All,listOfCovariancesOfH,sceneDataJointlyNormalised,eta);
    [listOfEstimatedH, diagnostic] = optimise_eta_peichen_algorithm(...
                            A,b,vs,ws,X_All,listOfCovariancesOfH,...
                                        sceneDataJointlyNormalised,eta);
  
    [dummy, numberOfPlanes] = size(listOfEstimatedH);
    
    for i = 1:numberOfPlanes
        listOfEstimatedH{i} =  Tp \ listOfEstimatedH{i} * T;
    end

end




% function [listOfT listOfTp] = normalise_individual_transform(sceneData)
%     numOfH = length(sceneData);
%     listOfT = cell(1,numOfH);
%     listOfTp = cell(1,numOfH);  
%     
%     for i = 1:numOfH
%         x = sceneData(i).ptsInViewOne;
%         xp = sceneData(i).ptsInViewTwo;
%         n = length(x);
%         x = [ x; ones( 1, n ) ];
%         xp = [ xp; ones( 1, n ) ];
%         T = normalize_transform( x );
%         Tp = normalize_transform( xp );
%         listOfT{i} = T;
%         listOfTp{i} = Tp;
%     end
% 
% end
% 
% function [sceneData listOfInitialH T Tp] = normalise_joint_transform(sceneData, listOfInitialH)
%     numOfH = length(sceneData);
% 
%     
%     % collect all the pts in view one into a matrix
%     % collect all the pts in view two into a matrix
%     x = sceneData(1).ptsInViewOne;
%     xp = sceneData(1).ptsInViewTwo;
%     for i = 2:numOfH
%       x = horzcat(x,sceneData(i).ptsInViewOne);
%       xp = horzcat(xp,sceneData(i).ptsInViewTwo);
%     end
%     
%     nn = length(x);
%     x = [ x; ones( 1, nn ) ];
%     xp = [ xp; ones( 1, nn ) ];
%     T = normalize_transform( x );
%     Tp = normalize_transform( xp );
%     
%     % convert the points and homographys in each plane into the new
%     % joint coordinate system    
%     for i = 1:numOfH
%         %listOfInitialH{i} = Tp \ listOfInitialH{i} * T;
%         listOfInitialH{i} = Tp * listOfInitialH{i} / T;
%         x = sceneData(i).ptsInViewOne;
%         xp = sceneData(i).ptsInViewTwo;
%         n = length(x);
%         x = [ x; ones( 1, n ) ];
%         xp = [ xp; ones( 1, n ) ];
%         
%         xx =  T * x;
%         xpp = Tp * xp;
%         
%         xx = xx ./ repmat( xx(3,:), 3, 1 );
%         xpp = xpp ./ repmat( xpp(3,:), 3, 1 );
%         sceneData(i).ptsInViewOne = xx(1:2,:);
%         sceneData(i).ptsInViewTwo = xpp(1:2,:);       
%     end
% 
% end

% function T = normalize_transform( x )
% 
%   % Transform taking x's centroid to the origin
% 
%   Ttrans = [ 1 0 -mean( x(1,:) ) ; 0 1 -mean( x(2,:) ) ; 0 0 1 ];
% 
%   % Calculate appropriate scaling factor
% 
%   x = Ttrans * x;
%   lengths = sqrt( sum( x(1:2,:).^2 ));
%   s = sqrt(2) / mean(lengths);
% 
%   % Transform scaling x to an average length of sqrt(2)
% 
%   Tscale = [ s 0 0 ; 0 s 0 ; 0 0 1 ];
% 
%   % Compose the transforms
% 
%   T = Tscale * Ttrans;
% end

function [listOfValidOptimisedH, diagnostic] = optimise_eta_peichen_algorithm(A,b,vs,ws,X_All,listOfCovariancesOfH,sceneDataJointlyNormalised,etap)
   
    % data structure used to keep track of diagnostic information relating
    % to the performance of the optimization method (e.g. running time etc)
    diagnostic.iterations = 0;
    diagnostic.residualHistory = [];
    diagnostic.runningTime = 0;
tic

maxIterations = 10000;
numberOfPlanes = length(vs);

% Prepare the covariance matrices
% This will be a block-diagonal matrix of covariances
Q = [];
listOfBi = cell(1,3);
for i = 1:numberOfPlanes
    x_i = X_All(:,i);
    covariance_i = listOfCovariancesOfH{i};
    
    % symmetric projection matrix
    Pox_i = eye(9,9) - (x_i * x_i')/ (norm(x_i,2)^2) ;
    
    Lambda = Pox_i * covariance_i * Pox_i;
    [U, S, V] = svd(Lambda);
    for n = 1:8
        S(n,n) = 1/(S(n,n))  ;
        
    end
    S(9,9) = 0;
    
    % Enforce rank-8 constraint on covariance
    Ai =  V * S *  U';
    
    % Take square root of covariance
    [U, S, V] = svd(Ai);
    for n = 1:9
        S(n,n) = sqrt(S(n,n))  ;
    end
      
    Bi = U*S*V'; 
    listOfBi{i} = Bi;
    Q = blkdiag(Q,Bi);
end


%--------------------------------------------
% H_All = [];
% vtemp = reshape(vs,3,numberOfPlanes);
% listOfTempH = cell(0,4);
% for c = 1:4
%     vn_k = vtemp(:,c);
%     
%     Hest_i = (ws(c)*A + b * vn_k');
%     %Hest_i = Hest_i / norm(Hest_i,'fro');
%     vech = vec(Hest_i);
%     %vech = vech / norm(vech)^2;
%     listOfTempH{c} = vech;
%     %H_All = [H_All Hest_i];
%     H_All = [H_All vech];
% end
% 
% cost_temp = 0;
% for c = 1:4
%     x_i = X_All(:,c);
%     x_i = x_i / norm(x_i)^2;
%     theta_i = listOfTempH{c};
%     1;
%     Bi = listOfBi{c};
%     %cost_temp = cost_temp + (theta_i - x_i)'*Bi*(theta_i - x_i);
%     %(theta_i )'*Bi*(theta_i )
%     %cost_temp = cost_temp + (theta_i )'*Bi*(theta_i );
%     % Corresponds with AML
%     %cost_temp = cost_temp + (theta_i )'*Bi'*Bi*(theta_i );
%     cost_temp = cost_temp + (theta_i - x_i)'*(Bi'*Bi)*(theta_i - x_i);
% end
% 1;
% 
% XN_All = X_All;
% HN_All = H_All;
% for c = 1:4
% XN_All(:,c) = XN_All(:,c) / norm(XN_All(:,c),2);
% HN_All(:,c) = HN_All(:,c) / norm(HN_All(:,c),2);
% end
% 
% cost(1) = (vec(HN_All) - vec(XN_All))' * (Q'*Q) * ( vec(HN_All) - vec(XN_All));
% eta2 = [vec(A)' vec(b)' vtemp(:,1)' ws(1) vtemp(:,2)' ws(2) vtemp(:,3)' ws(3) vtemp(:,4)' ws(4)]';
% 
% myAMLCost = AMLCost();
% [matrixOfCostOfEachPlane jacobianMatrix]= myAMLCost.cost_function_vector_aml(eta2,vec(X_All),listOfCovariancesOfH);
% cost_aml(1) = sum(sum(matrixOfCostOfEachPlane.^2));
% %------------------------------
1;

b_old = b;
a_old = vec(A);
w_old = ws;
v_old = vec(vs)

for iter = 1:3
    [w_new, v_new] = apply_post_processing_separate(a_old,b_old,w_old,v_old,sceneDataJointlyNormalised,numberOfPlanes);
    v_new = vec(v_new);
    [a_new, b_new] = apply_post_processing_variant(a_old,b_old,w_new,v_new,sceneDataJointlyNormalised,numberOfPlanes);
    a_old = a_new;
    b_old = b_new;
    w_old = w_new;
    v_old = v_new;
end


%cost_true = zeros(1,5);
%cost_pei = zeros(1,5);
% keeps track of the residuals for each iterations of the optimisation
residualHistory = [];
%tic    
counter = 1;
for i = 1:maxIterations
  if (i > 1)
    prevresnorm = resnorm;
  end
if (counter == 1)
    v_new = pinv(Q*( kron(eye(3*numberOfPlanes),b_old) )) * (Q * (vec(X_All)- kron(w_old,a_old)) );
    diff_v = sum((v_new-v_old).^2);
    v_old = v_new;
    1;
end

if (counter == 2)
   b_new = pinv(Q*( kron(v_old,eye(3)) ))  * (Q * (vec(X_All)- kron(w_old,a_old)) ) ;
   diff_b =  sum((b_new-b_old).^2);
   b_old = b_new;
   1;
end

if (counter == 4)    
    a_new = pinv(Q*( kron(w_old,eye(9)) ))  * (Q * (vec(X_All)- kron(eye(3*numberOfPlanes),b_old)*v_old) );
    diff_a = sum((a_new - a_old).^2);
    a_old = a_new;
    1;
end

if (counter == 3)
    w_new = pinv(Q*( kron(eye(numberOfPlanes),a_old) ))  * (Q * (vec(X_All)- kron(eye(3*numberOfPlanes),b_old)*v_old) );
    diff_w = sum((w_new - w_old).^2);
    w_old = w_new;
    1;
end
1;
%difference = sum((v_new-v_old).^2) + sum((b_new-b_old).^2) 
if (i > 4)
%     sum((v_new-v_old).^2)
%     sum((b_new-b_old).^2)
%     sum((a_new - a_old).^2)
%     sum((w_new - w_old).^2)
    diff_v;
    diff_b;
    diff_a;
    diff_w;
end
1;

%%%%%%%%%%%%%%%%%%%%%%% SANITY CHECK!!! %%%%%%%%%%%%%%%
H_All = [];
vtemp = reshape(v_new,3,numberOfPlanes);
%listOfTempH = cell(0,4);
for c = 1:numberOfPlanes
    vn_k = vtemp(:,c);
    A = reshape(a_new,3,3);
    Hest_i = (w_new(c)*A + b_new * vn_k');
    %Hest_i = Hest_i / norm(Hest_i,'fro');
    vech = vec(Hest_i);
    %vech = vech / norm(vech)^2;
    %listOfTempH{c} = vech;
    %H_All = [H_All Hest_i];
    H_All = [H_All vech];
end

resnorm = (vec(H_All) - vec(X_All))' * (Q'*Q) * ( vec(H_All) - vec(X_All));
residualHistory = [residualHistory; resnorm];

if (i > 1)
   if ( abs(resnorm - prevresnorm) < 1e-10)
    break;
   end
end
% 
% % cost_temp = 0;
% % for c = 1:4
% %     x_i = X_All(:,c);
% %     x_i = x_i / norm(x_i)^2;
% %     theta_i = listOfTempH{c};
% %     1;
% %     Bi = listOfBi{c};
% %     cost_temp = cost_temp + (theta_i - x_i)'*(Bi'*Bi)*(theta_i - x_i);
% % end
% 1;
% 
% XN_All = X_All;
% HN_All = H_All;
% for c = 1:4
% XN_All(:,c) = XN_All(:,c) / norm(XN_All(:,c),2);
% HN_All(:,c) = HN_All(:,c) / norm(HN_All(:,c),2);
% end

% cost_true(i) = (vec(HN_All) - vec(XN_All))' * (Q'*Q) * ( vec(HN_All) - vec(XN_All));
%cost_pei(i) = (vec(H_All) - vec(X_All))' * (Q'*Q) * ( vec(H_All) - vec(X_All));
% %eta2 = [vec(A)' vec(b)' vtemp(:,1)' ws(1) vtemp(:,2)' ws(2) vtemp(:,3)' ws(3) vtemp(:,4)' ws(4)]';
% eta2 = [a_new' b_new' vtemp(:,1)' w_new(1) vtemp(:,2)' w_new(2) vtemp(:,3)' w_new(3) vtemp(:,4)' w_new(4)]';
% myAMLCost = AMLCost();
% [matrixOfCostOfEachPlane jacobianMatrix]= myAMLCost.cost_function_vector_aml(eta2,vec(X_All),listOfCovariancesOfH);
% cost_aml(i) = sum(sum(matrixOfCostOfEachPlane.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%% SANITY CHECK %%%%%%%%%%%%%%%



counter = counter +1;
if (counter == 5)
    counter = 1;
end
% 1;

end


% diff_v
% diff_b
% diff_a
% diff_w
1;



1;
%[w_new vs] = apply_post_processing(a_new,b_new,w_new,v_new,sceneDataJointlyNormalised,numberOfPlanes);

%[w_new vs] = apply_post_processing_separate;

% works on bs and as
%[a_new b_new] = apply_post_processing_variant(a_new,b_new,w_new,v_new,sceneDataJointlyNormalised,numberOfPlanes);

%-----------------------------------------------------------
% Construct the estimated homographies from the parameters
listOfValidOptimisedH = cell(1,3);
A = reshape(a_new,3,3);

% USED TO BE before posthoc
vs = reshape(v_new,3,numberOfPlanes);

%%%%%%%%%%%%

b = b_new;
ws = w_new;
1;
% Compute homorgaphies
for k = 1:numberOfPlanes  
    vn_k = vs(:,k);     
    listOfValidOptimisedH{k} = (ws(k)*A + b * vn_k');
end

runningTime = toc; % end timer  

% store the diagnostic information 
diagnostic.iterations = i;
diagnostic.residualHistory = residualHistory;
diagnostic.runningTime = runningTime;

1;

end












function listOfValidOptimisedH = optimise_eta_peichen_algorithm_projective(A,b,vs,ws,X_All,listOfCovariancesOfH,sceneDataJointlyNormalised,etap)
maxIterations = 200;
%numberOfPlanes = length(vs);
[dummy, numberOfPlanes] = size(vs);

% Prepare the covariance matrices
% This will be a block-diagonal matrix of covariances
Q = [];
listOfBi = cell(1,2);
for i = 1:numberOfPlanes
    x_i = X_All(:,i);
    covariance_i = listOfCovariancesOfH{i};
    
    % symmetric projection matrix
    Pox_i = eye(9,9) - (x_i * x_i')/ (norm(x_i,2)^2) ;
    
    Lambda = Pox_i * covariance_i * Pox_i;
    [U, S, V] = svd(Lambda);
    for n = 1:8
        S(n,n) = 1/(S(n,n))  ;
        
    end
    S(9,9) = 0;
    
    % Enforce rank-8 constraint on covariance
    Ai =  V * S *  U';
    
    % Take square root of covariance
    [U, S, V] = svd(Ai);
    for n = 1:9
        S(n,n) = sqrt(S(n,n))  ;
    end
      
    Bi = U*S*V'; 
    listOfBi{i} = Bi;
    Q = blkdiag(Q,Bi);
end

1;

b_old = b;
a_old = vec(A);
w_old = ws;
v_old = vec(vs)

for iter = 1:3
    [w_new, v_new] = apply_post_processing_separate(a_old,b_old,w_old,v_old,sceneDataJointlyNormalised,numberOfPlanes);
    v_new = vec(v_new);
    [a_new, b_new] = apply_post_processing_variant(a_old,b_old,w_new,v_new,sceneDataJointlyNormalised,numberOfPlanes);
    a_old = a_new;
    b_old = b_new;
    w_old = w_new;
    v_old = v_new;
end


cost_true = zeros(1,5);
cost_pei = zeros(1,5);
counter = 1;
for i = 1:maxIterations
if (counter == 1)    
    v_new = pls_better(Q , kron(eye(3*numberOfPlanes),b_old) , kron(w_old,a_old),v_old);
    diff_v = sum((v_new-v_old).^2);
    v_old = v_new;
    1;
end

if (counter == 2)   
   b_new = pls_better(Q, kron(v_old,eye(3)),kron(w_old,a_old),b_old );
   diff_b =  sum((b_new-b_old).^2);
   b_old = b_new;
   1;
end

if (counter == 4)      
    a_new = pls_better(Q,kron(w_old,eye(9)),kron(eye(3*numberOfPlanes),b_old)*v_old,a_old);
    diff_a = sum((a_new - a_old).^2);
    a_old = a_new;
    1;
end

if (counter == 3)    
    w_new = pls_better(Q,kron(eye(numberOfPlanes),a_old),kron(eye(3*numberOfPlanes),b_old)*v_old,w_old );
    diff_w = sum((w_new - w_old).^2);
    w_old = w_new;
    1;
end
1;
%difference = sum((v_new-v_old).^2) + sum((b_new-b_old).^2) 
if (i > 4)
%     sum((v_new-v_old).^2)
%     sum((b_new-b_old).^2)
%     sum((a_new - a_old).^2)
%     sum((w_new - w_old).^2)
    diff_v;
    diff_b;
    diff_a;
    diff_w;
end
1;

%%%%%%%%%%%%%%%%%%%%%%% SANITY CHECK!!! %%%%%%%%%%%%%%%
%  H_All = [];
%  vtemp = reshape(v_new,3,numberOfPlanes);
% listOfTempH = cell(0,4);
% for c = 1:4
%     vn_k = vtemp(:,c);
%     A = reshape(a_new,3,3);
%     Hest_i = (w_new(c)*A + b_new * vn_k');
%     %Hest_i = Hest_i / norm(Hest_i,'fro');
%     vech = vec(Hest_i);
%     %vech = vech / norm(vech)^2;
%     listOfTempH{c} = vech;
%     %H_All = [H_All Hest_i];
%     H_All = [H_All vech];
% end
% 
% 1;
% 
% XN_All = X_All;
% HN_All = H_All;
% for c = 1:4
% XN_All(:,c) = XN_All(:,c) / norm(XN_All(:,c),2);
% HN_All(:,c) = HN_All(:,c) / norm(HN_All(:,c),2);
% end
% 
% cost_true(i) = (vec(HN_All) - vec(XN_All))' * (Q'*Q) * ( vec(HN_All) - vec(XN_All));
% cost_pei(i) = (vec(H_All) - vec(X_All))' * (Q'*Q) * ( vec(H_All) - vec(X_All));
%  eta2 = [a_new' b_new' vtemp(:,1)' w_new(1) vtemp(:,2)' w_new(2) vtemp(:,3)' w_new(3) vtemp(:,4)' w_new(4)]';
%  myAMLCost = AMLCost();
% [matrixOfCostOfEachPlane jacobianMatrix]= myAMLCost.cost_function_vector_aml(eta2,vec(X_All),listOfCovariancesOfH);
%  cost_aml(i) = sum(sum(matrixOfCostOfEachPlane.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%% SANITY CHECK %%%%%%%%%%%%%%%



counter = counter +1;
if (counter == 5)
    counter = 1;
end
% 1;
end

diff_v
diff_b
diff_a
diff_w
1;



1;


%-----------------------------------------------------------
% Construct the estimated homographies from the parameters
listOfValidOptimisedH = cell(1,2);
A = reshape(a_new,3,3);

% USED TO BE before posthoc
vs = reshape(v_new,3,numberOfPlanes);

%%%%%%%%%%%%

b = b_new;
ws = w_new;
1;
% Compute homorgaphies
for k = 1:numberOfPlanes  
    vn_k = vs(:,k);     
    listOfValidOptimisedH{k} = (ws(k)*A + b * vn_k');
end

1;

end

function y = pls(Z, A, b, y)
[k k2] = size(A);

J = norm(Z*(A*y+b))^2/ norm(A*y+b)^2;
y = pinv(A'*(Z'*Z-J*eye(k))*A)*A'*(J*eye(k) - Z'*Z)*b;


end

function y = pls_better(Z, A, b, y)
K = Z*[A b];
L = [A b];
[U,V,X,C,S] = gsvd(K,L);
X = inv(X)';
p = X(:,1);
y = p(1:end-1) / p(end);

end

% function y = pls_better(Z, A, b, y)
% K = Z*[A b];
% L = [A b];
% [V D] = eig(K'*K,L'*L);
% d = diag(D);
% [val ind] = sort(d);
% 
% p = V(:,ind(1));
% 
% 
% y = p(1:end-1) / p(end);
% 
% end


% Best one so far
function [w_fix, vs_fix] = apply_post_processing_separate(a_new,b_new,w_new,v_new,sceneDataJointlyNormalised,numberOfPlanes)
% Apply pei-chens post-processing

Ulist = cell(0,numberOfPlanes);
numberOfPointsArray = zeros(1,numberOfPlanes);
for k = 1:numberOfPlanes
    x = sceneDataJointlyNormalised(k).ptsInViewOne;
    xp = sceneDataJointlyNormalised(k).ptsInViewTwo;
    n = length(x);
    numberOfPointsArray(k) = n;
    U = [];
    % Construct block diagonal data carrier matrices
    for i=1:n
        % the x-coordinate in the first image
        u = x(1,i);
        % the y-coordinate in the first image
        v = x(2,i);
        % the x-coordinate in the second image
        u2 = xp(1,i);
        % the y-coordinate in the second image
        v2 = xp(2,i);
        
        %         % Constructor Data Carrier matrix
        %         ux1 = [ 0    0    0  -u  -v -1  u*v2   v*v2   v2 ]';
        %         ux2 = [ u    v    1   0   0  0  -u*u2 -v*u2  -u2 ]';
        %         ux3 = [-u*v2 -v*v2 -v2 u*u2 v*u2 u2  0    0      0 ]';
        %         % 9 x 3 matrix
        %         Ui = [ux1 ux2 ux3];
        
        m = [u v 1]';
        m2 = [u2 v2 1]';
        M =[0 -m2(3) m2(2); m2(3) 0 -m2(1); -m2(2) m2(1) 0];
        Ui = kron(-m,M);
        U = vertcat(U,Ui');
    end
    Ulist{k} = U;
    
end

S = [kron(eye(3),b_new) a_new];
T =  [reshape(v_new,3,numberOfPlanes); w_new'];
1;

Tnew = zeros(4,numberOfPlanes);
for i = 1:numberOfPlanes    
[W, D, V] = svd(Ulist{i}*S);
Tnew(:,i) = V(:,end);
1;
end
1;
%%%%%%%%%%%
w_fix = (Tnew(4,:))';
vs_fix = Tnew(1:3,:);
end



function [w_fix, vs_fix] = apply_post_processing(a_new,b_new,w_new,v_new,sceneDataJointlyNormalised,numberOfPlanes)
% Apply pei-chens post-processing
U = [];
numberOfPointsArray = zeros(1,numberOfPlanes);
for i = 1:numberOfPlanes
    x = sceneDataJointlyNormalised(i).ptsInViewOne;
    xp = sceneDataJointlyNormalised(i).ptsInViewTwo;
    n = length(x);
    numberOfPointsArray(i) = n;
    % Construct block diagonal data carrier matrices
    for i=1:n
        % the x-coordinate in the first image
        u = x(1,i);
        % the y-coordinate in the first image
        v = x(2,i);
        % the x-coordinate in the second image
        u2 = xp(1,i);
        % the y-coordinate in the second image
        v2 = xp(2,i);
        
        %         % Constructor Data Carrier matrix
        %         ux1 = [ 0    0    0  -u  -v -1  u*v2   v*v2   v2 ]';
        %         ux2 = [ u    v    1   0   0  0  -u*u2 -v*u2  -u2 ]';
        %         ux3 = [-u*v2 -v*v2 -v2 u*u2 v*u2 u2  0    0      0 ]';
        %         % 9 x 3 matrix
        %         Ui = [ux1 ux2 ux3];
        
        m = [u v 1]';
        m2 = [u2 v2 1]';
        M =[0 -m2(3) m2(2); m2(3) 0 -m2(1); -m2(2) m2(1) 0];
        Ui = kron(-m,M);
        U = blkdiag(U,Ui);
    end
    
end

fArray = zeros(sum(numberOfPointsArray),numberOfPlanes);
F = [];
for i = 1:numberOfPlanes
    f = zeros(sum(numberOfPointsArray),1);
    onesVec = ones(numberOfPointsArray(i),1);
    1;
    if i == 1
        f(1:numberOfPointsArray(i)) = onesVec;
        fArray(:,i)= f;
        F = [F kron(f,eye(9))];
        1;
    else
        starting = sum(numberOfPointsArray(1:i-1))+1;
        ending = starting + numberOfPointsArray(i);
        f(starting:ending-1) = onesVec;
        fArray(:,i)= f;
        F = [F kron(f,eye(9))];
        1;
    end
    
    1;
end
1;
S = [kron(eye(3),b_new) a_new];
T =  [reshape(v_new,3,numberOfPlanes); w_new'];
1;
[U, S, V] = svd(U'*F*kron(eye(numberOfPlanes),S));
1;
%%%%%%%%%%%
Tnew = V(:,end-3)+V(:,end-2) + V(:,end-1) + V(:,end)
Tnew = reshape(Tnew,4,4)

w_fix = (Tnew(4,:));
vs_fix = Tnew(1:3,:);
end

function [a_fix, b_fix] = apply_post_processing_variant(a_new,b_new,w_new,v_new,sceneDataJointlyNormalised,numberOfPlanes)
% Apply pei-chens post-processing
U = [];
numberOfPointsArray = zeros(1,numberOfPlanes);
for i = 1:numberOfPlanes
    x = sceneDataJointlyNormalised(i).ptsInViewOne;
    xp = sceneDataJointlyNormalised(i).ptsInViewTwo;
    n = length(x);
    numberOfPointsArray(i) = n;
    % Construct block diagonal data carrier matrices
    for i=1:n
        % the x-coordinate in the first image
        u = x(1,i);
        % the y-coordinate in the first image
        v = x(2,i);
        % the x-coordinate in the second image
        u2 = xp(1,i);
        % the y-coordinate in the second image
        v2 = xp(2,i);
        
        %         % Constructor Data Carrier matrix
        %         ux1 = [ 0    0    0  -u  -v -1  u*v2   v*v2   v2 ]';
        %         ux2 = [ u    v    1   0   0  0  -u*u2 -v*u2  -u2 ]';
        %         ux3 = [-u*v2 -v*v2 -v2 u*u2 v*u2 u2  0    0      0 ]';
        %         % 9 x 3 matrix
        %         Ui = [ux1 ux2 ux3];
        
        m = [u v 1]';
        m2 = [u2 v2 1]';
        M =[0 -m2(3) m2(2); m2(3) 0 -m2(1); -m2(2) m2(1) 0];
        Ui = kron(-m,M);
        U = blkdiag(U,Ui);
    end
    
end

fArray = zeros(sum(numberOfPointsArray),numberOfPlanes);
F = [];
for i = 1:numberOfPlanes
    f = zeros(sum(numberOfPointsArray),1);
    onesVec = ones(numberOfPointsArray(i),1);
    1;
    if i == 1
        f(1:numberOfPointsArray(i)) = onesVec;
        fArray(:,i)= f;
        F = [F kron(f,eye(9))];
        1;
    else
        starting = sum(numberOfPointsArray(1:i-1))+1;
        ending = starting + numberOfPointsArray(i);
        f(starting:ending-1) = onesVec;
        fArray(:,i)= f;
        F = [F kron(f,eye(9))];
        1;
    end
    
    1;
end
1;
S = [kron(eye(3),b_new) a_new];
T =  [reshape(v_new,3,numberOfPlanes); w_new'];
1;
z = U'*F*kron(T',eye(9))*vec(S);
1;
K1 = [kron(vec(eye(3)),eye(3)); zeros(9,3)];
K2 = [zeros(27,9) ; eye(9)];
K = [K1 K2];
k = [b_new; a_new];
1;
[U, S, V] = svd(U'*F*kron(T',eye(9))*K);
1;
%%%%%%%%%%%
k = V(:,end);
b_fix = k(1:3);
a_fix = k(4:end);
1;
end
