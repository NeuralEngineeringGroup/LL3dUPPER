function Dataout=Near_NaN_Euclidian(Data0,Number_K_near,graph)
%%

[Np,Ndim,Ns] = size(Data0);
Data = reshape(Data0,[Np*Ndim,Ns]);
Dataout = Data;
   
incomplete_samples = find(any(~isfinite(Data),1)); % return nan index of data
for ii=1:length(incomplete_samples)
    Y=Data(:,incomplete_samples(ii));
    missing_y = ~isfinite(Y); % return the nan index of vector sample 
    
    nonmissing_samples = all(isfinite(Data(missing_y,:)),1);
    Data_reduced = Data(:,nonmissing_samples);
    
    N_reduced=size(Data_reduced,2);
    Y0=Y*ones(1,N_reduced);

    dif2=(Data_reduced-Y0).^2;
    [~,indsort] = sort(mean(dif2,"omitnan"),'ascend');
    indsort_asc_KNN = indsort(1:Number_K_near);
    X = Data_reduced(:,indsort_asc_KNN);
    T = mean(X,2,"omitnan");

    Y(missing_y) = T(missing_y);
    Dataout(:,incomplete_samples(ii)) = Y;
    %ii
    
end
Dataout = reshape(Dataout,Np,Ndim,Ns);

if graph
    figure
    subplot(1,2,1)
    imagesc(isnan(Data))
    %xlim([1 1000])
    title('Original Data')
    ylabel('3D-Poses')
    xlabel('Number of Sample')
    %
    subplot(1,2,2)
    imagesc(isnan(Dataout))
    %xlim([1 1000])
    title('Recounstructed Data Base on KNN')
    ylabel('3D-Poses')
    xlabel('Number of Sample')
end
%%
end
