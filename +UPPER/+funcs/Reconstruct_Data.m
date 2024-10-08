function [Data_reconstruct_3D]=Reconstruct_Data(RawData3D_full,Data_3D_KNN,Threshold_Outliers,mean_pose_3D,mean_pose_ppca,Cov_pPCA)
%%

[Np, Framedim,Nsample]=size(RawData3D_full);
%Data_3D_KNN=RawData3D_full;
is_outlier = false(Np, Framedim,Nsample);
for n = 1:Nsample
    is_outlier(:,:,n) = UPPER.funcs.detect_outliers(squeeze(Data_3D_KNN(:,:,n)), mean_pose_3D, Cov_pPCA, Threshold_Outliers);  
end
Data_3D_KNN(is_outlier==1) = NaN;
% Data_2D_KNN = reshape(Data_3D_KNN,Np*Framedim,Nsample);
% Outlier_percent_fram = (length(find(sum(isnan(Data_2D_KNN))))/(Nsample))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%re-align data without outliers (very important!!!) 
Data_3D_alignment_WO = UPPER.funcs.Alignment(Data_3D_KNN, mean_pose_3D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%reconstruction
Data_3D_alignment_WO_reshape = reshape(Data_3D_alignment_WO, Np*Framedim,Nsample);
Data_reconstruct = UPPER.funcs.theoretical_estimate_correction(Data_3D_alignment_WO_reshape,mean_pose_ppca,Cov_pPCA);
Data_reconstruct_3D = reshape(Data_reconstruct, Np, Framedim, Nsample);
Data_reconstruct_3D_align = UPPER.funcs.Alignment(Data_reconstruct_3D, mean_pose_3D);
%back to original place
RawData3D=RawData3D_full;
Reconstructed_data_final = zeros(size(RawData3D_full),"Like",RawData3D_full); % Preallocation is a guess, might be wrong size
for jj = 1:length(RawData3D_full)
	kk = find((sum(isnan(RawData3D(:,:,jj)),2))>1);
	Data_reconstruct=Data_reconstruct_3D_align(:,:,jj);
	Dra_Raw=RawData3D(:,:,jj);
	if isempty(kk)==1
		[~, ~, Trans1] = procrustes(Dra_Raw, Data_reconstruct, 'Scaling', false,'Reflection',false);
		Y_sample2 = Trans1.b*Data_reconstruct_3D_align(:,:,jj)*Trans1.T + repmat(Trans1.c(1,:),Np,1);
		Reconstructed_data_final(:,:,jj) = Y_sample2;
    	
	%         reconst=Alignment_ric(Data_reconstruct,Dra_Raw);
	%         Reconstructed_data_final(:,:,j)=reconst;
	else
		Data_reconstruct(kk,:,:)=[];
		Dra_Raw(kk,:,:)=[]; 
		[~, ~, Trans] = procrustes(Dra_Raw, Data_reconstruct, 'Scaling', false,'Reflection',false);
		Y_sample = Trans.b*Data_reconstruct_3D_align(:,:,jj)*Trans.T + repmat(Trans.c(1,:),Np,1);
		Reconstructed_data_final(:,:,jj) = Y_sample;
	end
end
Data_reconstruct_3D=Reconstructed_data_final;    

%%
end
