function [points, new_points] = pointsRefine(points, marked)
    % sparse grid magic to find neighbors of marked points
    [old_points, d] = size(points.t); marked = reshape(marked, 1, []);
    t = inf(2 * d * numel(marked), d); level = t; new_points = 0;
    for j = marked
        if points.used(j)
            warning("skipping already used point")
        else
            for k = 1:d
                dt = points.dt(k);
                % add a point to the "left" and "right" in direction k
                for dir = [-1 1]
                    t_new = points.t(j, :); level_new = points.level(j, :);
                    level_new(:, k) = level_new(:, k) + 1;
                    t_new(k) = t_new(k) + dir * dt * .5^level_new(k);
                    if t_new(k) >= points.t(1, k) - dt && t_new(k) <= points.t(1, k) + dt
                        % check if hierarchy hash already exists
                        h_found = find(all(points.hierarchy_hash == level_new, 2));
                        t_already_exists = false; % whether point exists already!
                        if h_found
                            % look into right bin and compare to see if point exists
                            idx_hash = points.hierarchy{h_found};
                            idx_hash_old = idx_hash(idx_hash <= old_points);
                            if any(all(points.t(idx_hash_old, :) == t_new, 2))
                                t_already_exists = true;
                            else
                                idx_hash_new = idx_hash(idx_hash > old_points) - old_points;
                                if any(all(t_new == t(idx_hash_new, :), 2))
                                    t_already_exists = true;
                                end
                            end
                        end
                        if ~t_already_exists
                            new_points = new_points + 1;
                            t(new_points, :) = t_new;
                            level(new_points, :) = level_new;
                            if h_found % append to old hash
                                points.hierarchy{h_found}(end + 1) = old_points + new_points;
                            else % create new hash
                                points.hierarchy_hash(end + 1, :) = level_new;
                                points.hierarchy{end + 1} = old_points + new_points;
                            end
                        end
                    end
                end
            end
            points.used(j) = 1;
        end
    end
    points.t = [points.t; t(1 : new_points, :)];
    points.level = [points.level; level(1 : new_points, :)];
    points.used = [points.used zeros(1, new_points)];
end
