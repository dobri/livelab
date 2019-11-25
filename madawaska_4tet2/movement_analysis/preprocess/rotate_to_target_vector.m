function DATA = rotate_to_target_vector(DATA,target_vector,dims,bodies,reference_markers)

clf
% Do this separately for each body.
for b = 1:numel(bodies)
    % Find the indices of the needed markers.
    marker_label{1} = [bodies{b} reference_markers{1}];
    marker_label{2} = [bodies{b} reference_markers{2}];
    marker_index = zeros(1,2);
    for m = 1:numel(DATA.col_names)
        for marker = 1:2
            if strcmp(DATA.col_names{m},marker_label{marker})
                marker_index(marker) = m;
            end
        end
    end
    
    % Skip ahead and don't do anything if you didn't find all markers.
    if any(marker_index==0)
        %keyboard
        continue
    end
    
    % Extract the 3D values of these markers.
    X = DATA.X(:,dims,marker_index);
    
    % Find the local vector A->B frame by frame.
    local_vector = zeros(size(X,1),numel(dims));
    for k = 1:size(X,1)
        local_vector(k,:) = X(k,:,2)-X(k,:,1);
    end
    
    % Then the angle from the target vector.
    angles_vec = zeros(size(X,1),1);
    Rz = zeros(3,3);
    for k = 1:size(X,1)
        if local_vector(k,2)>0
            angles_vec(k) = -acos(dot(local_vector(k,:),target_vector)/...
                (norm(local_vector(k,:))*norm(target_vector)));
        else
            angles_vec(k) = acos(dot(local_vector(k,:),target_vector)/...
                (norm(local_vector(k,:))*norm(target_vector)));
        end
        %Rz(:,:,k) = [cos(-angles_vec(k)) -sin(-angles_vec(k)) 0; sin(-angles_vec(k)) cos(-angles_vec(k)) 0;0 0 1];
    end
    % Define the rotation matrix. Assume we're rotation about Z.
    Rz = [cos(nanmean(-angles_vec)) -sin(nanmean(-angles_vec)) 0; sin(nanmean(-angles_vec)) cos(nanmean(-angles_vec)) 0;0 0 1];
    % Check that the angle makes sense.
    % Very roughly the angles should be somewhere around 60, 30, -60, -80 
    % for violin1, 2, viola, and cello, respectively.
    subplot(numel(bodies),1,b)
    plot(angles_vec/2/pi*360)
    
    % Now get ready to rotate. 
    % First, find the indices of the markers that belong to the given body.
    markers_ind_per_body = false(1,numel(DATA.col_names));
    for m = 1:numel(DATA.col_names)
        if strcmp(DATA.col_names{m}(1:numel(bodies{b})),bodies{b})
            markers_ind_per_body(m) = true;
        end
    end
    X = DATA.X(:,:,markers_ind_per_body);
    
    % Zero-center the body.
    center_vec = nanmean(X,3);
    X = X-center_vec;
    % Rotate frame by frame.
    for k = 1:size(X,1)
        for m = 1:size(X,3)
            %X(k,:,m) = (Rz(:,:,k)*X(k,:,m)')';
            X(k,:,m) = (Rz*X(k,:,m)')';
        end
    end
    % Un-zero-center the body.
    X = X+center_vec;
    DATA.X(:,:,markers_ind_per_body) = X;
end

