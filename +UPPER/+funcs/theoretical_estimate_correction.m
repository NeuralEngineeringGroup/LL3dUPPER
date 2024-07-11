function Data_Estimate_L=theoretical_estimate_correction(data,mean_data,Cov_pPCA)
Data_Estimate_L=data;
for ii=1:size(data,2)
    Y_sample=data(:,ii);
    if sum(isnan(Y_sample))>0
        Index_nan = ~isfinite(Y_sample);
        Y_sample_E2=Y_sample;
        Y_sample_E2(Index_nan)=[];
        mean_E2=mean_data;
        mean_E1=mean_data(Index_nan,:);
        mean_E2(Index_nan)=[];
        %
        Cov_pPCA_E22=Cov_pPCA(~Index_nan,~Index_nan);
        %
        Cov_pPCA_E12=Cov_pPCA(Index_nan,isfinite(Y_sample));
        Mu_Estimate= mean_E1+Cov_pPCA_E12*pinv(Cov_pPCA_E22)*(Y_sample_E2-mean_E2);
        Y_sample(Index_nan,:)=Mu_Estimate;
        Data_Estimate_L(:,ii)=Y_sample;
    end
end   