function DATA = pca_all_trials_in_DATA(DATA)

% We want these markers on these bodies.
bodies_labels = {'violin1','violin2','viola','cello'};
headMark = {'cellohat0','cellohat1','cellohat2','cellohat3','violahat0','violahat1','violahat2','violahat3',...
    'violin1hat0','violin1hat1','violin1hat2','violin1hat3','violin2hat0','violin2hat1','violin2hat2','violin2hat3'};

% Find and average the four head markers per body. Then PCA, and keep some information.
for b = 1:numel(bodies_labels)
    for tr = 1:numel(DATA)
        markerInd = [];
        for m = 1:numel(headMark)
            if ~isempty(regexp(headMark{m},[bodies_labels{b} '(\w*)'],'once'))
                markerInd = [markerInd find(strcmp(DATA{tr}.col_names, headMark{m}))];
            end
        end
        X = nanmean(DATA{tr}.X_clean(:,:,markerInd),3);
        [P,Xrotated,~,~,r2] = pca(X);
        DATA{tr}.Xpcs(:,:,b) = Xrotated;
        DATA{tr}.r2(:,b) = r2';
        [~,ind] = max(P,[],2);
        DATA{tr}.dim_mapped_mostly_to(:,b) = ind';
        fprintf('Trial %1.0f, body: %8s, R2 per PC: %5.2f%%, %5.2f%%, %5.2f%%. Which dims: %1.f, %1.f, %1.f.\n',tr,bodies_labels{b},r2,ind)
    end
    %{
            Inspect P to see which dimensions of X go where after the transformation.
            Column n of P tells you how much the dimensions of the original
            data are mapped to dimension n after PCA. So if you have three
            dimensions of input and the first column has a much higher value in its
            second value, then this means that the second dimension of the
            original data is mapped heavily onto the first dimension after
            transformation. Kapish?
        
            For example, the mapping
            P = [0.0823    0.9905   -0.1098
                 0.9893   -0.0944   -0.1109
                 0.1202    0.0995    0.9878];
            means that the first dimension of the PCAed data is made almost
            enitrely of the second original dimension. Respectively,
            the second of the first. And the third of the third.
            
            To get a feel how it works, follow the steps here while also
            looking at the results. First we plot the output of PCA. Then we
            "manually" reconstruct again the first dimension of this output
            and plot it on top. You see we have exactly recreated the 1st d
            of the transformed data by doing matrix multiplication with the
            first column of P.
            plot(Xrotated);hold on;plot((P(:,1)'*(X-mean(X))')','o');hold off
    %}
end
