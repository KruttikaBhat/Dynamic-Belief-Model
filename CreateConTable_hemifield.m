function [contable, contable_lc, contable_rc] = CreateConTable_hemifield(data, all_angles)

num_angles = length(all_angles);
contable = zeros(3,3,num_angles);
contable_lc = zeros(3,3,num_angles);
contable_rc = zeros(3,3,num_angles);

% 2: valid; 3: invalid; 1: no change.
validity_column = [2,3,1];

% 1:overall (leftcue and rightcue combined); 2:leftcue; 3:rightcue 
for cue_type = 1:3
    for block = 1: size(data,3)
        for angle = 1: num_angles
            ctable = zeros(3,3);
            
            if cue_type == 1
                idx = find(data(:,1,block) ~= 0 & data(:,9,block) ~= 5 & data(:,4,block) == all_angles(angle) );
                idx_nc = find(data(:,1,block) ~= 0 & data(:,9,block) ~= 5 & data(:,2,block) == 1 );
                
            elseif cue_type == 2
                idx = find(data(:,1,block) ~= 0 & data(:,9,block) ~= 5 & data(:,4,block) == all_angles(angle)...
                    & data(:,3,block) == -1);
                idx_nc = find(data(:,1,block) ~= 0 & data(:,9,block) ~= 5 & data(:,2,block) == 1 ...
                    & data(:,3,block) == -1);
                
            else
                idx = find(data(:,1,block) ~= 0 & data(:,9,block) ~= 5 & data(:,4,block) == all_angles(angle)...
                    & data(:,3,block) == 1);
                idx_nc = find(data(:,1,block) ~= 0 & data(:,9,block) ~= 5 & data(:,2,block) == 1 ...
                    & data(:,3,block) == 1);
            end
            
            for row = 1:2
                val = validity_column(row);
                ctable(row,1) = length( find(data(idx,2,block) == val & data(idx,9,block) == data(idx,3,block)) );
                ctable(row,2) = length( find(data(idx,2,block) == val & data(idx,9,block) == -data(idx,3,block)) );
                ctable(row,3) = length( find(data(idx,2,block) == val & data(idx,9,block) == 0) );
            end
            
            ctable(3,1) = length( find(data(idx_nc,9,block) == data(idx_nc,3,block)) );
            ctable(3,2) = length( find(data(idx_nc,9,block) == -data(idx_nc,3,block)) );
            ctable(3,3) = length( find(data(idx_nc,9,block) == 0) );
            
            if cue_type == 1
                contable(:,:,angle) = contable(:,:,angle) + ctable;
            elseif cue_type == 2
                contable_lc(:,:,angle) = contable_lc(:,:,angle) + ctable;
            else
                contable_rc(:,:,angle) = contable_rc(:,:,angle) + ctable;
            end
            
        end
    end
end
