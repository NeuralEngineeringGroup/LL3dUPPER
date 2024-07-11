function [data_filled] = knn_imputation(data,num_neighbors,plot_results)
%%
arguments
	data			(:,:,:)	{mustBeFloat}
	num_neighbors	(1,1)	double
	plot_results	(1,1)	logical		= false
end
%%
data_filled = data;
fill_idx = find(any(~isfinite(data),[1 2]));
for ii = fill_idx(:).'
	d = data(:,:,ii);
	nonmissing_samples = all(isfinite(d) | isfinite(data), [1 2]);
	sample_dist = mean((data - d).^2,[1 2],"omitnan");
	sample_quantile = num_neighbors / sum(nonmissing_samples);
	dist_threshold = quantile(sample_dist(nonmissing_samples),sample_quantile);
	matching_samples = nonmissing_samples(:) & sample_dist(:)<=dist_threshold;
	Template = mean(data(:,:,matching_samples),3,"omitnan");

	d(~isfinite(d)) = Template(~isfinite(d));
	data_filled(:,:,ii) = d;
end

if plot_results
    figure
    subplot(1,2,1)
    imagesc(isnan(reshape(data,[],size(data,3))))
    %xlim([1 1000])
    title('Original Data')
    ylabel('Coordinates')
    xlabel('Number of Sample')
    %
    subplot(1,2,2)
    imagesc(isnan(reshape(data_filled,[],size(data,3))))
    %xlim([1 1000])
    title('Recounstructed Data Base on KNN')
    ylabel('Coordinates')
    xlabel('Number of Sample')
end
%%
end
