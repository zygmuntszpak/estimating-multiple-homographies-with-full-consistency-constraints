%   Function: visualise_symmetric_transfer_errors_for_fixed_noise
%
%   Computes the average homography estimation error for various homography estimation methods.
%   An average error is first computed for each scene. The average is taken by averaging over
%   all homographies/planes in a scene for each scene. Once all the average errors in a scene
%   are computed, a further average is taken over all scenes, resulting in a single number quantifying 
%   the performance of homography estimation methods. 
%
%   Parameters:
%
%      listOfErrorsForEachMethodAndScene     -   a cell array containing the symmetric transfer errors
%											     for each scene, each homography estimation method and 
%											     each homography/plane.
%
%   Returns: 
% 
%   A cell array quantifying the average performance of each homography estimation method.
%
%   See Also:
%
%  compute_symmetric_transfer_error_for_all_methods
%  visualise_mean_symmetric_transfer_errors_for_all_noise
%  
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

function listOfMeanErrorsForEachMethod = visualise_symmetric_transfer_errors_for_fixed_noise( listOfErrorsForEachMethodAndScene )
    
    [numberOfScenes numberOfEstimationMethods numberOfPlanes] = size(listOfErrorsForEachMethodAndScene);
    rmse_dlt = zeros(numberOfScenes,numberOfPlanes);
    rmse_fns = zeros(numberOfScenes,numberOfPlanes);
    rmse_aml_fns = zeros(numberOfScenes,numberOfPlanes);
    rmse_aml_dlt = zeros(numberOfScenes,numberOfPlanes);
    rmse_ba_dlt_1 = zeros(numberOfScenes,numberOfPlanes);
    rmse_ba_fns_1 = zeros(numberOfScenes,numberOfPlanes);
    rmse_chen_fns = zeros(numberOfScenes,numberOfPlanes);
    rmse_chen_dlt = zeros(numberOfScenes,numberOfPlanes);
    rmse_aml_smps_fns = zeros(numberOfScenes,numberOfPlanes);

    for i=1:numberOfScenes
        for j=1:numberOfPlanes
 
        % Since each point is represented by an x and y coordinate we
        % divide the total number of entries in the vector by 2 
        numberOfPoints = length(listOfErrorsForEachMethodAndScene{i,NDLT_CONST,j}.errorsForEachPlane) / 2;
        
        rmse_dlt(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,NDLT_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));
        1;
        rmse_fns(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,FNS_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));
        rmse_aml_fns(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,AML_FNS_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));
        rmse_aml_dlt(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,AML_DLT_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));      
        rmse_ba_dlt_1(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,BA_DLT_1_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));       
        rmse_ba_fns_1(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,BA_FNS_1_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));   
        rmse_chen_fns(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,CHEN_FNS_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));
        rmse_chen_dlt(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,CHEN_DLT_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));    
        rmse_aml_smps_fns(i,j) = sqrt(sum(listOfErrorsForEachMethodAndScene{i,AML_SMPS_FNS_CONST,j}.errorsForEachPlane)/ (numberOfPlanes*numberOfPoints));   

        end
    end
    
   

    
    fprintf('\n ERRORS FOR DLT: \n')
    
    rmse_dlt;
    
    % averaged across all scenes
    mean_rmse_dlt = mean(rmse_dlt)
    
    % averaged across all planes
    mean_mean_rmse_dlt = mean(mean_rmse_dlt)
    
    fprintf('\n ERRORS FOR FNS: \n')
    
    rmse_fns;
    
    % averaged across all scenes
    mean_rmse_fns = mean(rmse_fns)
    
    % averaged across all planes
    mean_mean_rmse_fns = mean(mean_rmse_fns)
    
    
     fprintf('\n ERRORS FOR AML-FNS: \n')
    
    rmse_aml_fns;
    
    % averaged across all scenes
    mean_rmse_aml_fns = mean(rmse_aml_fns)
    
    % averaged across all planes
    mean_mean_rmse_aml_fns = mean(mean_rmse_aml_fns)
    
   
    
    fprintf('\n ERRORS FOR AML-DLT: \n')
    
    rmse_aml_dlt
    
    % averaged across all scenes
    mean_rmse_aml_dlt = mean(rmse_aml_dlt)
    
    % averaged across all planes
    mean_mean_rmse_aml_dlt = mean(mean_rmse_aml_dlt)
    
    
    fprintf('\n ERRORS FOR BA-DLT-1: \n')
    
    rmse_ba_dlt_1;
    
    % averaged across all scenes
    mean_rmse_ba_dlt_1 = mean(rmse_ba_dlt_1)
    
    % averaged across all planes
    mean_mean_rmse_ba_dlt_1 = mean(mean_rmse_ba_dlt_1)
    
  
    fprintf('\n ERRORS FOR BA-FNS-1: \n')
    
    rmse_ba_fns_1;
    
    % averaged across all scenes
    mean_rmse_ba_fns_1 = mean(rmse_ba_fns_1)
    
    % averaged across all planes
    mean_mean_rmse_ba_fns_1 = mean(mean_rmse_ba_fns_1)

    fprintf('\n ERRORS FOR CHEN-FNS: \n')
    
    rmse_chen_fns
    
    % averaged across all scenes
    mean_rmse_chen_fns = mean(rmse_chen_fns)
    
    % averaged across all planes
    mean_mean_rmse_chen_fns = mean(mean_rmse_chen_fns)
    
    
    fprintf('\n ERRORS FOR CHEN-DLT: \n')
    
    rmse_chen_dlt
    
    % averaged across all scenes
    mean_rmse_chen_dlt = mean(rmse_chen_dlt)
    
    % averaged across all planes
    mean_mean_rmse_chen_dlt = mean(mean_rmse_chen_dlt)
    
    fprintf('\n ERRORS FOR AML-SMPS-FNS: \n')
    
    rmse_aml_smps_fns
    
    % averaged across all scenes
    mean_rmse_aml_smps_fns = mean(rmse_aml_smps_fns)
    
    % averaged across all planes
    mean_mean_rmse_aml_smps_fns = mean(mean_rmse_aml_smps_fns)
	
	% Uncomment the code below to display boxplots comparing the performance of the various estimation methods.
    
%     %--------------------------------
%     figure
%     boxplot([rmse_dlt(:,1) rmse_fns(:,1) rmse_aml_dlt(:,1) rmse_aml_fns(:,1)  rmse_chen_dlt(:,1) rmse_chen_fns(:,1)   rmse_ba_dlt_1(:,1) rmse_ba_fns_1(:,1) ],'notch','on','boxstyle','outline','plotstyle','traditional','labelorientation','inline','labels',{'DLT';'FNS'; 'AML-DLT';'AML-FNS'; 'WALS-DLT';'WALS-FNS'; 'BA-DLT';'BA-FNS'})
%     
%     figure
%     boxplot([rmse_dlt(:,2) rmse_fns(:,2) rmse_aml_dlt(:,2) rmse_aml_fns(:,2)  rmse_chen_dlt(:,2) rmse_chen_fns(:,2)   rmse_ba_dlt_1(:,2) rmse_ba_fns_1(:,2)],'notch','on','boxstyle','outline','plotstyle','traditional','labelorientation','inline','labels',{'DLT';'FNS'; 'AML-DLT';'AML-FNS'; 'WALS-DLT';'WALS-FNS'; 'BA-DLT';'BA-FNS'})
%     
%     figure
%     boxplot([rmse_dlt(:,3) rmse_fns(:,3) rmse_aml_dlt(:,3) rmse_aml_fns(:,3)  rmse_chen_dlt(:,3) rmse_chen_fns(:,3)   rmse_ba_dlt_1(:,3) rmse_ba_fns_1(:,3)],'notch','on','boxstyle','outline','plotstyle','traditional','labelorientation','inline','labels',{'DLT';'FNS'; 'AML-DLT';'AML-FNS'; 'WALS-DLT';'WALS-FNS'; 'BA-DLT';'BA-FNS'})
%     
%     figure
%     boxplot([rmse_dlt(:,4) rmse_fns(:,4) rmse_aml_dlt(:,4) rmse_aml_fns(:,4)  rmse_chen_dlt(:,4) rmse_chen_fns(:,4)   rmse_ba_dlt_1(:,4) rmse_ba_fns_1(:,4)],'notch','on','boxstyle','outline','plotstyle','traditional','labelorientation','inline','labels',{'DLT';'FNS'; 'AML-DLT';'AML-FNS'; 'WALS-DLT';'WALS-FNS'; 'BA-DLT';'BA-FNS'})
%     %--------------------------------
   

    listOfMeanErrorsForEachMethod = cell(1,9);
    listOfMeanErrorsForEachMethod{NDLT_CONST} = mean_mean_rmse_dlt;
    listOfMeanErrorsForEachMethod{FNS_CONST} =  mean_mean_rmse_fns; 
    listOfMeanErrorsForEachMethod{AML_FNS_CONST} = mean_mean_rmse_aml_fns;        
    listOfMeanErrorsForEachMethod{AML_DLT_CONST} = mean_mean_rmse_aml_dlt;
    listOfMeanErrorsForEachMethod{BA_DLT_1_CONST} = mean_mean_rmse_ba_dlt_1;
    listOfMeanErrorsForEachMethod{BA_FNS_1_CONST} = mean_mean_rmse_ba_fns_1;
    listOfMeanErrorsForEachMethod{CHEN_FNS_CONST} = mean_mean_rmse_chen_fns;
    listOfMeanErrorsForEachMethod{CHEN_DLT_CONST} = mean_mean_rmse_chen_dlt;
    listOfMeanErrorsForEachMethod{AML_SMPS_FNS_CONST} = mean_mean_rmse_aml_smps_fns;
    
end

