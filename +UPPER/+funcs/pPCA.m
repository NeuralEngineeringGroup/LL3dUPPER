function [mean_pose_vec, NumDimcut, Cov_pPCA, eignValues, eignVectors]=pPCA(data,Threshold_Eigen,plot_results,opts)
%ref
%https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.934.5867&rep=rep1&type=pdf
%%
arguments
	data				(:,:,:)	{mustBeFloat}
	Threshold_Eigen		(1,1)	double	= 0.95
	plot_results		(1,1)	logical	= false
	opts.cov_opts		(1,1)	struct = struct(Method="ogk",NumOGKIterations=2)
end
%%

[Nparts,Ndims,N]=size(data);
data_vec = reshape(data,[Nparts*3,N]);
NP = Nparts*Ndims;

mean_pose_vec = mean(data_vec,2,"omitnan"); % transpose Posev to calculate mean from all frame for x,y,z all 11 poses and then again transpose.

cov_opts = namedargs2cell(opts.cov_opts);
Cov_Data0 = robustcov((data_vec - mean_pose_vec).', cov_opts{:});

%%%% PCA
[eignVectors,eignValues] = eig(Cov_Data0,'vector');
[eignValues,indsort] = sort(eignValues,'descend');
eignVectors = eignVectors(:,indsort);

%construct the optimal hyperplane with the error projection
error_project = zeros(length(eignValues),1);
for k=1:length(eignValues)
    error_project(k)=sum(eignValues(1:k))/sum(eignValues);
end

r0=find(error_project>Threshold_Eigen); % obtain cut off eigenvector
NumDimcut=r0(1);
sigma2=mean(eignValues(r0(1)+1:end)); % averaging from remaining eigenvalues


%%% pPCA
diagonal_vector=[eignValues(1:r0(1))',sigma2*ones(1,NP-r0(1))]; %% defining Sigma^2 
Cov_pPCA=eignVectors*(diag(diagonal_vector))*eignVectors';

if plot_results
    %
    figure 
    pcolor(Cov_Data0-Cov_pPCA)
    colorbar
    title('Cov-Original - Cov-pPCA')
    %
    figure
    plot(error_project,'LineWidth',3)
    xlabel('#Eigenvalues')
    ylabel('Variance Expl')
    ylim([0. 1.01])
    xticks(1:1:length(eignValues))
    xline(r0(1),'--g','LineWidth',2)
end

%%
end