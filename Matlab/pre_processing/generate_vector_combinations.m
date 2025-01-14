%GENERATE_VECTOR_COMBINATIONS Generates a matrix (v_combo) with all
%   possible combinations of the simulation parameters entered as a vector. In
%   addition the associated column index for each vector entry is stored in
%   v_index.
%
% Copyright (c) 2021 Hannah Dentzien, Ove Hagge Ellhöft
% Copyright (c) 2021 Jens Geisler

function [DLC_cell, v_combo, v_index, g_index] = generate_vector_combinations(DLC_cell, row_xls, col_start, config)
%% Identify all vectors in row & save all possible combinations

v_combo = [1];        % intialize vector combination matrix to use allcombos function
v_index = [];         % storage for column indices of vectors in DLC_cell
g_index = 0;         % variation group number
struct_id = struct(); % storage for vectors with identifiers

% allow to define Vektor with respect to reference wind conditions
v_i= str2double(config.CutinWind);
v_r= str2double(config.RatedWind);
v_o= str2double(config.CutoutWind);

% loop over each "non-basic" column
for col_xls = col_start:size(DLC_cell,2)
    if strcmp(DLC_cell{row_xls,col_xls}, 'v_i')
        DLC_cell{row_xls,col_xls}= v_i;
    elseif strcmp(DLC_cell{row_xls,col_xls}, 'v_r')
        DLC_cell{row_xls,col_xls}= v_r;
    elseif strcmp(DLC_cell{row_xls,col_xls}, 'v_o')
        DLC_cell{row_xls,col_xls}= v_o;
    end
    
    % create all combinations of vectors WITHOUT indentifiers
    try  
        e = eval(DLC_cell{row_xls,col_xls}); % read vector from char
        if isvector(e) && numel(e)>1      % check if element is vector with more than one element
            v_combo = allcombos(v_combo,e); % combination of vectors
            v_index = [v_index, col_xls]; % save row number of rows with vector
            g_index(end+1)= g_index(end) + 1;
        end
    catch
        % process vectors only; try next column
    end

    % special treatment for uni wind conditions: a list of conditions
    % separated by colons
    if strcmp(DLC_cell{1, col_xls}, '{uni-wind-param}')
        n_iec= sum(DLC_cell{row_xls, col_xls}==':')+1;
        if n_iec>1
            v_combo = allcombos(v_combo, 1:n_iec); % combination of vectors
            v_index = [v_index, col_xls]; % save row number of rows with vector            
            g_index(end+1)= g_index(end) + 1;
        end
    end

    % sort all vectors WITH identifiers by used identifier
    if ischar(DLC_cell{row_xls,col_xls}) && contains(DLC_cell{row_xls,col_xls},'<') % check if identifier is used
        split = strsplit(DLC_cell{row_xls,col_xls},'>'); % split element into identifier and vector
        ident_name = strip(split{1},'left','<');         % save name of used identifier

        if isfield(struct_id,ident_name)   % check if identifier has already been used

            if length(struct_id.(ident_name))-1 == length(eval(split{2})) % vectors with identifiers have to be of same legth
                struct_id.(ident_name)(end+1, :) = [col_xls, eval(split{2})]; % add vector to exsisting identifier; save column of vector as first element in each row of struct
            else 
                error('DLC %s: Dimensions of vectors with identifier <%s> are not consistent.\n', DLC_cell{row_xls,1},ident_name);
            end

        else
            struct_id.(ident_name) = [col_xls, eval(split{2})]; % create new field with identifier
        end
    end

end

% create combinations of ALL vectors 
ids = fieldnames(struct_id);   % check if struct contains identifiers 
for idx = 1:length(ids)
   v_combo = allcombos(v_combo, struct_id.(ids{idx})(:,2:end));  % combine identifier-vectors with other vectors
   groups= (struct_id.(ids{idx})(:,1))';
   v_index = [v_index, groups];            % save columns of identifier-vectors
   g_index = [g_index, (g_index(end)+1)*ones(size(groups))];
end

v_combo = v_combo(2:end,:);    % erase initial value of v_combo
g_index = g_index(2:end);
