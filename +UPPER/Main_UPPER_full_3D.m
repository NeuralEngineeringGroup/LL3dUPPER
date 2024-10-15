% Main Body code for the estimation of Stastical Shape Model(SSM)

min_p = 0.8;
trial_xyz = load("trial_xyz.mat").trial_xyz;

% while sum(cellfun(@height, XYZ.xyz))>10000
% 	XYZ = XYZ(1:2:end,:);
% end

XYZ_P = cat(1,trial_xyz.p{:});
XYZ = cat(1,trial_xyz.xyz{:});
XYZ = XYZ + 0./all(XYZ_P>=min_p,3);
clear XYZ_P

XYZ = permute(XYZ, [2 3 1]);
XYZ([2 13 14 15],:,:) = [];

%%
mdl = UPPER.UPPER3dModel();
idx = find(sum(~all(isfinite(XYZ),2),1)<=3);
idx = sort(idx(randperm(numel(idx),5e4)));
mdl.data = XYZ(:,:,idx);

%%%
% Estimation_Model function receives the 3D data and return the Data which all missing data filled up and  mean
% estimated by Ransac, Mean of pPCA, Covariance of pPCA and Eigenvalues
% and eigenpose
% SSM = struct;
% [SSM.knn, SSM.ransac_mean, SSM.ppca_mean, SSM.ppca_cov, SSM.eigenvalues, SSM.eigenvectors] ...
% 	=UPPER.funcs.Estimation_Model(XYZ_sub,thresh=0.8);


tic
mdl = mdl.fit(thresh=0.8);
toc
%%
%the input is the 3D Rawdata (it is suggested use the Data_3D_KNN which has no missing values) and the output is 3D reconstructed data which
%backed to original place in the arena. 

% [Recounstructed_Data_full]=UPPER.funcs.Reconstruct_Data(...
% 	RawData3D_full, ...
% 	SSM.knn, ...
% 	0.99, ...
% 	SSM.ransac_mean, ...
% 	SSM.ppca_mean, ...
% 	SSM.ppca_cov)
tic
XYZr = mdl.reconstruct(XYZ(:,:,1:37:end),outlier=0.99);
toc
%%
% save('Recounstructed_Data_full')


