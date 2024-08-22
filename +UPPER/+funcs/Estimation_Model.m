function [Data_3D_KNN, mean_pose_3D, mean_pose_ppca, Cov_pPCA, eignValues, eignVectors]=Estimation_Model(RawData3D,opts)
arguments (Input)
	RawData3D					(:,3,:)
	opts.thresholdEigenValue	(1,1)	double	= 0.8
	opts.plotResults			(1,1)	logical	= false
end
%this function receives data including missing data and outliers. estimate
%the robust covariance matrix and return the covariance matrix and eigenpose
%and eigenvalues
%%%%
% [Np,Framedim,Ns] = size(RawData3D); % Unused
mean_pose_3D = UPPER.funcs.estimate_mean_pose_RANSAC(RawData3D, false);
%%%%%%%
%alignment
Data_3D_align = UPPER.funcs.Alignment(RawData3D, mean_pose_3D);
%%%%%%%
%Filling data with KNN=5
K = 5;
Data_3D_KNN = UPPER.funcs.knn_impute_points(Data_3D_align, K, opts.plotResults);
% Data_3D_KNN_P = Data_3D_KNN; % Unused
% Data_KNN_reshape = reshape(Data_3D_KNN, Np*Framedim,Ns); % Unused
%%%%%%%
%PPCA
[mean_pose_ppca, ~, Cov_pPCA, eignValues, eignVectors] = UPPER.funcs.pPCA(Data_3D_KNN,opts.thresholdEigenValue,true);
%mean_pose_3D_ppca = reshape(mean_pose_ppca,[Np,Framedim]);

end
