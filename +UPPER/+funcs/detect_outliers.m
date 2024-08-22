function[is_outlier] = detect_outliers(X,mu,C,outlier_threshold,opts)
%%
arguments
	X
	mu
	C
	outlier_threshold			(1,1)	double	= 0.99
	opts.threshold_multiplier	(1,1)	double
end
%%

[Nbp,Ndim] = size(X);
if ~isfield(opts,"threshold_multiplier")
	if Ndim==2,	opts.threshold_multiplier = 25;
	else,		opts.threshold_multiplier = 3;
	end
end

inlier  = find(~isnan(X(:,1)));
Ninlier = numel(inlier);
stop_search = false;
while ~stop_search
	%generate full list of body point subsets
	list_bp = nchoosek(inlier, Ninlier);
	Nlist = size(list_bp,1);
	%initialise square mahalanobis distance
	dist = zeros(1,Nlist);
	for n = 1:Nlist
		%get subset of body points, mean and covariance
		Xsub = X(list_bp(n,:),:); 
		mean_pose_sub = mu(list_bp(n,:),:);
		if Ndim==3
			indCsub = [list_bp(n,:) list_bp(n,:)+Nbp list_bp(n,:)+2*Nbp]; % 3D version
		else
			indCsub = [list_bp(n,:) list_bp(n,:)+Nbp]; % 2D version
		end
		Csub_inv = C(indCsub, indCsub)^-1;

		%re-align pose for the subset
		[~, Xsub] = procrustes(mean_pose_sub, Xsub, 'Scaling', false, 'Reflection', false);
		%calculate squared mahalanobis distance
		dist(n) = (Xsub(:)-mean_pose_sub(:))'*Csub_inv*(Xsub(:)-mean_pose_sub(:));
	end
	%find minimum distance
	[minC, indmin] = min(dist);
	%label body points associated with minC as inliers
	inlier = list_bp(indmin,:);
	
	%compare with threshold
	TH = opts.threshold_multiplier*chi2inv(outlier_threshold,Ndim*Ninlier);
	if minC<TH
		stop_search = true;
	elseif Ninlier<0.5*Nbp
		inlier = [];
		stop_search = true;
	else
		Ninlier = Ninlier-1;
	end 
end
is_outlier = true(Nbp,Ndim); 
is_outlier(inlier,:) = false;
%%
end