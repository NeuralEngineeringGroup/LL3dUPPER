function [RawData]=Alignment(RawData,Ref_Pose)
[Np, Ndim, Num_samples]=size(RawData);
for n = 1:Num_samples
    Y_sample = squeeze(RawData(:,:,n));
    indnum = find(~isnan(Y_sample(:,1)));
    if numel(indnum) > max(0.5*Np, Ndim)
        [~, ~, Trans] = procrustes(Ref_Pose(indnum,:), Y_sample(indnum,:), 'Scaling', false,'Reflection',false);
        Y_sample = Trans.b*Y_sample*Trans.T + repmat(Trans.c(1,:),Np,1);
        RawData(:,:,n) = Y_sample;
    else
        RawData(:,:,n) = NaN;
    end
end
