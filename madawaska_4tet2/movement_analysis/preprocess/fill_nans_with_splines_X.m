function [x_filled] = fill_nans_with_splines_X(D)
    x_filled = D;
    x_filled(1,:)=nanmean(D);
    x_filled(end,:)=nanmean(D);
    for d = 1:3
        % Fill small nan gaps by linearly interpolating. Be careful with larger gaps!
        for k=2:size(x_filled,1)-1
            if isnan(x_filled(k,d))
                k_end = min([size(D,1) k-1+find(~isnan(D(k:end,d)),1,'first')]);
                x_filled(k:k_end-1,d) = linspace(x_filled(k-1,d),x_filled(k_end,d),k_end-k);
            end
        end
    end
end