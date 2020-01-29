function x = fill_nans_by_lin_interp(x)

for d=1:size(x,2)
    % Fill small nan gaps by linearly interpolating. Be careful with larger gaps!
    x([1 end],d)=0;
    for k=2:size(x,1)-1
        if isnan(x(k,d))
            k_end = min([size(x,1) k-1+find(~isnan(x(k:end,d)),1,'first')]);
            x(k:k_end-1,d) = linspace(x(k-1,d),x(k_end,d),k_end-k);
        end
    end
end