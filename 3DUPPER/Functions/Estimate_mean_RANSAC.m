function mean_sample_3D=Estimate_mean_RANSAC(data,plot_results,opts)
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
New_ref = zeros(Nparts,Ndims,opts.iter,typ);
for ii=1:opts.iter
    %%% Select a random pose as estimate for mean
    Rand_ind_ref = randsample(index_Non_Nan,1);
    mean_sample_3D = data(:,:,Rand_ind_ref);
	
    Rand_ind_sample=randsample(N,sample_size);
    dataSample=data(:,:,Rand_ind_sample);
    %%%% Align
    Sample_Align_3D=Alignment(dataSample,mean_sample_3D);
	
	%%% Squared Euclidean distance
	Dist2 = mean((Sample_Align_3D - mean_sample_3D).^2,[1 2]);
    %%% Count neighbours
    NNeigh(ii) = sum(Dist2<dist_threshold);
    New_ref(:,:,ii) = mean(Sample_Align_3D(:,:,ind),3,"omitnan");
end
[~,best_candidate]=max(NNeigh);  %<---- Aghileh Change
mean_sample_3D=New_ref(:,:,best_candidate);

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