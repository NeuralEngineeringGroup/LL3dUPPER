function [mean_pose] = estimate_mean_pose_RANSAC(data,plot_results,opts)
%Riccardo Suggestion for alignment (Riccardo modified the 200421 code in
%Github)
arguments
	data			(:,:,:)	{mustBeFloat}
	plot_results	(1,1)	logical		= false	

	opts.iter				(1,1)	double		= 100
	opts.sample_ratio		(1,1)	double		= 0.25
	opts.num_rand_pairs		(1,1)	double		= 90-size(data,2)*20
	opts.neighbor_threshold	(1,1)	double		= 0.25
end
%% Get size & type of input data
[Nparts,Ndims,N]=size(data);
typ = class(data);
sample_size = round(N*opts.sample_ratio);
index_Non_Nan = find(all(isfinite(data),[1 2]));

%% Tune TH value based on 25th quantile of paiwise distances
dist = zeros(1,(opts.num_rand_pairs^2 - opts.num_rand_pairs)/2,typ);
dist_idx = 1;
Rand_ind = randsample(index_Non_Nan,opts.num_rand_pairs);
for n = 1:opts.num_rand_pairs-1
    temp = squeeze(data(:,:,Rand_ind(n)));
    for m = n+1:opts.num_rand_pairs
        temp1 = squeeze(data(:,:,Rand_ind(m)));
        [~,temp1] = procrustes(temp, temp1,'Scaling',false,'Reflection',false);
        dist(dist_idx) = mean((temp(:)-temp1(:)).^2,"omitnan");
		dist_idx = dist_idx+1;
    end
end
dist_threshold = quantile(dist,opts.neighbor_threshold);

%% RANSAC
NNeigh = zeros(opts.iter,1);
candidate_means = zeros(Nparts,Ndims,opts.iter,typ);
for ii=1:opts.iter
    candidate = data(:,:,randsample(index_Non_Nan,1)); % Select a random pose as estimate for mean
    sample = data(:,:,randsample(N,sample_size)); % Select random sample

    Sample_Align_3D = UPPER.funcs.Alignment(sample,candidate); % Align sample points to candidate
	
	Dist2 = mean((Sample_Align_3D - candidate).^2,[1 2]); % Squared Euclidean distance

	inliers = Dist2<dist_threshold; % Identify inliers (neighbors)

    NNeigh(ii) = sum(inliers); % Count neighbours
    candidate_means(:,:,ii) = mean(Sample_Align_3D(:,:,inliers),3,"omitnan"); % Compute new mean pose
end
[~,best_candidate] = max(NNeigh);  %<---- Aghileh Change
mean_pose = candidate_means(:,:,best_candidate);

%%% figure
if plot_results
    figure
    plot(NNeigh)
    title('Scoring of mean-pose with Close Poses');
    ylabel('Close Poses');
    xlabel('Iteration');
end
%%
end