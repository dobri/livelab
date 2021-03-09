function [SigmaA,SigmaV] = accelerations_after_smooth_and_threshold_cleaning(X,sf,t,plotting_flag)

% Extract the 3D values of these markers.
%X = ;

% We don't actually need all four head markers. %X = X(:,:,1);
% Try averaging across them to smooth the data a little.
X = nanmean(X,3);

dt = 1./sf;
% Also, verify that some markers do not disappear a lot, which
% would make the average jump up/down by a few mm.
% Some glitches in the raw 3D cause huge differences in v.
% It's easier to pre-smooth X before v, rather than to clean v.
V0 = [[0 0 0];diff(X)]./dt;
spikes = any(abs(V0)>1e3,2);
V = V0;
V(spikes,:)=nan;
for d=1:size(V0,2)
    V(:,d) = fill_nans_by_lin_interp(V(:,d));
    % Smooth by a third of a second. With the mov ave method this
    % kills everything above 3 Hz. With sgolay above 5 Hz.
    V(:,d) = smooth(V(:,d),round(sf/3),'sgolay');
end
SigmaV = dot(V,V,2).^.5;

A0 = [[0 0 0];diff(V)]./dt;
%A0 = [[0 0 0];diff(V0)]./dt; to see some of the discontinuities.
spikes = any(abs(A0)>10e3,2);
A = A0;
A(spikes,:)=nan;
for d=1:size(A,2)
    A(:,d) = fill_nans_by_lin_interp(A(:,d));
    % Smooth by a third of a second. With the mov ave method this
    % kills everything above 3 Hz. With sgolay above 5 Hz.
    A(:,d) = smooth(A(:,d),round(sf/3),'sgolay');
end
SigmaA = dot(A,A,2).^.5;

% Cut 3 secs in and out just in case.
SigmaA([1:300 end-299:end],:) = 0;

% Verify that there aren't spikes and nans remaining by accident.
% Did we do a good job cleaning and filtering without killing v?
% Compare v_temp (pure diff) w/ v after rem nans, lin interp, smooth.
if plotting_flag == 1
    figure(1)
    subplot(3,1,1)
    plot(t,X,'linewidth',2)
    subplot(3,1,2)
    plot(t,V0,'-','linewidth',2)
    hold on
    plot(t,V,'--','linewidth',2)
    hold off
    legend('V0x','V0y','V0z','Vx','Vy','Vz')
    subplot(3,1,3)
    plot(t,SigmaV,'-','linewidth',2)

    figure(2)
    subplot(2,1,1)
    plot(t,A0,'-','linewidth',2)
    hold on
    plot(t,A,'--','linewidth',2)
    hold off
    legend('A0x','A0y','A0z','Ax','Ay','Az')
    subplot(2,1,2)
    plot(t,SigmaA,'-','linewidth',2)
    hold on
    plot(t,spikes*max(SigmaA),'r-','linewidth',2)
    hold off
end