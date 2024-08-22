classdef UPPER3dModel
	properties
		data				(:,3,:)		single	= single.empty(0,3,0)
		bodyparts			(1,:)		string
		nearest_neighbors	(:,3,:)		single	{mustBeNonmissing}	= single.empty(0,3,0)
		mean_pose			(:,3)		single	{mustBeNonmissing}	= single.empty(0,3)

		PPCA_mean			(:,1)		single	{mustBeNonmissing}	= single.empty(0,1)
		PPCA_cov			(:,:)		single	{mustBeNonmissing}	= single.empty(0,0)

		eigenvalues			(:,1)		single	{mustBeNonmissing}	= single.empty(0,1)
		eigenvectors		(:,:)		single	{mustBeNonmissing}	= single.empty(0,0)
	end

	methods
		function obj = fit(obj,opts)
			arguments (Input)
				obj
				opts.thresholdEigenvalue
			end
			optcell = namedargs2cell(opts);

			[...
				obj.nearest_neighbors, ...
				obj.mean_pose, ...
				obj.PPCA_mean, ...
				obj.PPCA_cov, ...
				obj.eigenvalues, ...
				obj.eigenvectors ...
				] = UPPER.funcs.Estimation_Model(obj.data, optcell{:});
		end

		function reconstructed_data = reconstruct(obj, data, opts)
			arguments
				obj
				data					(:,3,:)	single
				opts.outlierThreshold	(1,1)	double	= 0.99
			end
			
			
			[Np, Framedim,Nsample]=size(data);
			%obj.nearest_neighbors=data;
			is_outlier = false(size(data));
			for n = 1:Nsample
    			is_outlier(:,:,n) = UPPER.funcs.detect_outliers(squeeze(data(:,:,n)), obj.mean_pose, obj.PPCA_cov, opts.outlierThreshold);  
			end
			data(is_outlier==1) = NaN;

			% Data_2D_KNN = reshape(obj.nearest_neighbors,Np*Framedim,Nsample);
			% Outlier_percent_fram = (length(find(sum(isnan(Data_2D_KNN))))/(Nsample))*100;
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%re-align data without outliers (very important!!!) 
			data_aligned = UPPER.funcs.Alignment(data, obj.mean_pose);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%reconstruction
			data_aligned = reshape(data_aligned, Np*Framedim,Nsample);
			data_aligned = UPPER.funcs.theoretical_estimate_correction(data_aligned,obj.PPCA_mean,obj.PPCA_cov);
			reconstruction_data = reshape(data_aligned, Np, Framedim, Nsample);
			reconstruction_data = UPPER.funcs.Alignment(reconstruction_data, obj.mean_pose);
			%back to original place
			
			reconstructed_data = NaN(size(data),"Like",data); % Preallocation is a guess, might be wrong size
			for jj = 1:size(data,3)
				kk = find((sum(isnan(data(:,:,jj)),2))>1);
				data_aligned=reconstruction_data(:,:,jj);
				Dra_Raw=data(:,:,jj);

				data_aligned(kk,:,:)=[];
				Dra_Raw(kk,:,:)=[];
				if isempty(data_aligned)
					continue
				end

				[~, ~, Trans] = procrustes(Dra_Raw, data_aligned, 'Scaling', false,'Reflection',false);
				reconstructed_data(:,:,jj) = Trans.b*reconstruction_data(:,:,jj)*Trans.T + repmat(Trans.c(1,:),Np,1);
			end
		end


	end
end