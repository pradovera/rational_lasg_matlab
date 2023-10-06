function [points] = pointsBuildHierarchy(points)
    % sparse grid magic to build hierarchical structure for sparse grid
    T = size(points.t, 1); points.hierarchy_hash = [points.level(1, :)]; points.hierarchy = {1};
    for j = 2:T
        % check if hierarchy hash already exists
        h_found = find(all(points.hierarchy_hash == points.level(j, :), 2));
        if h_found
            % append to old hash
            points.hierarchy{h_found}(end + 1) = j;
        else
            % create new hash
            points.hierarchy_hash(end + 1, :) = points.level(j, :);
            points.hierarchy{end + 1} = j;
        end
    end
end