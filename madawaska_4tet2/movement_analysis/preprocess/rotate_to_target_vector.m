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
    local_vectorAB_vec = X(:,:,2)-X(:,:,1);
    local_vectorAB = nanmean(local_vectorAB_vec);
    
    % Then the angle from the target vector.
    % Careful with trig functions in different quadrants!
    if local_vectorAB(2)>0
        the_angle = -acos(dot(local_vectorAB,target_vector)/...
            (norm(local_vectorAB)*norm(target_vector)));
    else
        the_angle = acos(dot(local_vectorAB,target_vector)/...
            (norm(local_vectorAB)*norm(target_vector)));
    end
    
    % For verification, also get a time series of angles.
    angles_vec = zeros(size(X,1),1);
    for k = 1:size(X,1)
        if local_vectorAB_vec(k,2)>0
            angles_vec(k) = -acos(dot(local_vectorAB_vec(k,:),target_vector)/...
                (norm(local_vectorAB_vec(k,:))*norm(target_vector)));
        else
            angles_vec(k) = acos(dot(local_vectorAB_vec(k,:),target_vector)/...
                (norm(local_vectorAB_vec(k,:))*norm(target_vector)));
        end
        % What to average? Ave the angles or vector and then compute one Rz
        % or many Rz and average them?
        %Rz(:,:,k) = [cos(-angles_vec(k)) -sin(-angles_vec(k)) 0; sin(-angles_vec(k)) cos(-angles_vec(k)) 0;0 0 1];
    end
    
    % Check that the angle makes sense. Very very roughly the angles t.s.
    % should be somewhere around 60, 30, -60, -80 for v1, v2, viola, and cello, respectively.
    subplot(numel(bodies),1,b)
    plot(angles_vec/2/pi*360)

    % Define the rotation matrix. Assume we're rotating about Z.
    Rz = [cos(-the_angle) -sin(-the_angle) 0;
        sin(-the_angle) cos(-the_angle) 0;
        0 0 1];
    
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
    
    % Plug the rotated markers back in the data struct.
    DATA.X(:,:,markers_ind_per_body) = X;
end

