path = 'Z:\EXPERIMENT-DATA\2022_degenerate_He3_and_He4_mixture\he3_he4_mixture';

out_folders = {{},{},{},{}};
data_folders = find_files(path,out_folders);

combined_struct = @(S,T) cell2struct(cellfun(@vert_or_horz_cat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);

%, types 'atom_num','flip','dead_laser_+_ai'

%% log labview
folder_list = data_folders{1};
shot_type_list = cell(1,length(folder_list));
all_variable = {'unsorted'};
for ii = 1:length(folder_list)
    current_dir = folder_list{ii};
    current_files = dir(current_dir);
    name_file = {current_files.name};

    log_mask = cellfun(@(x) contains(x,'log_'),name_file);
    if sum(log_mask)>0
        log_file = name_file(log_mask);
        log_file = log_file{1};
        log_table = readtable(fullfile(current_dir,log_file));
        variable_vec = log_table{:,4};
        if ~iscell(variable_vec(1))
            variable_vec=cellstr(num2str(variable_vec));
        end
        variable{ii} = unique(variable_vec);
        all_variable = [all_variable, setdiff(unique(variable_vec),all_variable).'];
        %     var_str = 'shot type';%'hold time';
        num_var = size(variable{ii});
        indx_list = log_table{:,1}.';%1:size(variable_vec,1);
        %     indx_list = log_table{:,1}.';%1:size(variable_vec,1);
        shot_list = cell(1,length(variable{ii}));
        for jj = 1:num_var
            shot_list{jj} = indx_list(strcmp(variable_vec,variable{ii}{jj}));%
        end
        shot_type_list{ii} = shot_list;
    else
        variable{ii} = 'unsorted';
        shot_type_list{ii} = {};
    end
end

%% import data
opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;
use_types = {'unsorted','  0.1','  0.2','  0.4','  0.8','  1.6','0.054','mix_shot'};
clear data
for ii = 1:length(folder_list)
    opts.import.dir = folder_list{ii};
    opts.import.cache_save_dir = fullfile(folder_list{ii}, 'cache', 'import\');

    if ~isempty(shot_type_list{ii})
        type_mask = contains(variable{ii},use_types);
        opts.import.shot_num = cell2mat(shot_type_list{ii}(type_mask)); %can select specific shots to import
        current_varaibles = variable{ii}(type_mask);
        current_shot_list = shot_type_list{ii}(type_mask);
    elseif isfield(opts.import,'shot_num')
        opts.import=rmfield(opts.import,'shot_num');
        current_varaibles = {'unsorted'};
        current_shot_list = {1:1e5};
    else
        current_varaibles = {'unsorted'};
        current_shot_list = {1:1e5};
    end
    fprintf('Setting up for:\n %s \n', folder_list{ii})
    
    if ii == 1
        [data, ~] = import_mcp_tdc_data(opts.import);
        data.dir = cell(size(data.counts_txy));
        data.dir(:)={folder_list{ii}};
        data.tag = cell(size(data.counts_txy));
        for jj = 1:length(current_varaibles)
        tga_mask = ismember(data.shot_num,current_shot_list{jj});
        data.tag(tga_mask)={current_varaibles{jj}};
        end
    else
        [data_temp, ~] = import_mcp_tdc_data(opts.import);
        data_temp.dir = cell(size(data_temp.counts_txy));
        data_temp.dir(:)={folder_list{ii}};
        data_temp.tag = cell(size(data_temp.counts_txy));
        for jj = 1:length(current_varaibles)
        tga_mask = ismember(data_temp.shot_num,current_shot_list{jj});
        data_temp.tag(tga_mask)={current_varaibles{jj}};
        end
        data = combined_struct(data,data_temp);
    end
end
%% adaptive time controlls
num_min = 0e3;
num_max = Inf;
he3_lim = [0,Inf];

he4_time = [-1,-1];
he3_time = [0,4];

shot_mask = data.num_counts>num_min & data.num_counts<num_max;

data_masked.num_counts = data.num_counts;
data_masked.counts_txy = data.counts_txy;
data_masked = struct_mask(data_masked,shot_mask);

%%
num_min = 50;

he4_times = {[0.437,0.46],[1.05,1.09]};
he3_times = {[0.41,0.437],[1.01,1.05]};

pal_time = [0.41,3.1];

he4_cen = [-3.7796,-4.16263].*1e-3;
he3_cen = [-3.827,-6.5065].*1e-3;

he3_lim = [10,25e3];%[500,1000];
he4_lim = [-1,15e3];

ax = 1; %which axis do we analyse [t,x,y or radial]

axis_labe = {'t','x','y','r'};
kk = 1;
clear he_dist
for ii = 1:length(data.counts_txy)
    if isempty(data.counts_txy{ii})
        continue
    end

    if size(masktxy_square(data.counts_txy{ii},[0.41 0.46; -0.03, 0.03; -0.03, 0.03]),1)>...
            size(masktxy_square(data.counts_txy{ii},[1.01 1.09; -0.03, 0.03; -0.03, 0.03]),1)
        he4_time = he4_times{1};
        he3_time = he3_times{1};
    else
        he4_time = he4_times{2};
        he3_time = he3_times{2};
    end
    
    lims_4 = [he4_time; -0.03, 0.03; -0.03, 0.03];
    he4_txy = masktxy_square(data.counts_txy{ii}, lims_4);
    lims_3 = [he3_time; -0.03, 0.03; -0.03, 0.03];
    he3_txy = masktxy_square(data.counts_txy{ii}, lims_3);

    t_he4 = linspace(he4_time(1),he4_time(2),5e3).';
    t_he3 = linspace(he3_time(1),he3_time(2),5e3).';

    space_bins = linspace(-30e-3,30e-3,5e3).';
    rad_bins = linspace(0.1e-3,30e-3,5e3).';

    % all_txy = ;

    if ax == 1
        bin_cen_4 = t_he4;
        bin_cen_3 = t_he3;
    elseif ax ~= 4
        bin_cen_4 = space_bins;
        bin_cen_3 = space_bins;
    else
        bin_cen_4 = rad_bins;
        bin_cen_3 = rad_bins;
    end

    %% run good shot checks
    is_shot_good = (size(he4_txy,1)+size(he3_txy,1))>num_min;
    he3_num_check = size(he3_txy(:,1),1)>he3_lim(1) && size(he3_txy(:,1),1)<he3_lim(2);
    he4_num_check = size(he4_txy(:,1),1)>he4_lim(1) && size(he4_txy(:,1),1)<he4_lim(2);
    is_shot_good= is_shot_good&he3_num_check&he4_num_check;
    shot_check(ii) = is_shot_good;
    if is_shot_good
        %% histogram in time
        sigma = 0.5e-4;
        if ax == 4
            count_hist_he4 = smooth_hist(sqrt((he4_txy(:,2)-he4_cen(1)).^2+(he4_txy(:,3)-he4_cen(2)).^2),'sigma',sigma,'edges',bin_cen_4);
            count_hist_he4.count_rate.smooth = count_hist_he4.count_rate.smooth./(2.*pi*count_hist_he4.bin.centers);
            count_hist_he3 = smooth_hist(sqrt((he3_txy(:,2)-he3_cen(1)).^2+(he3_txy(:,3)-he3_cen(2)).^2),'sigma',sigma,'edges',bin_cen_3);
            count_hist_he3.count_rate.smooth = count_hist_he3.count_rate.smooth./(2.*pi*count_hist_he3.bin.centers);
        else
            count_hist_he4 = smooth_hist(he4_txy(:,ax),'sigma',sigma,'edges',bin_cen_4);
            count_hist_he3 = smooth_hist(he3_txy(:,ax),'sigma',sigma,'edges',bin_cen_3);

            %take average over radial direction
            % sigmar = 0.8e-4;
            % count_hist_he4 = smooth_hist(sqrt((he4_txy(:,2)-he4_cen(1)).^2+(he4_txy(:,3)-he4_cen(2)).^2),'sigma',sigmar,'edges',rad_bins);
            % count_hist_he3 = smooth_hist(sqrt((he3_txy(:,2)-he3_cen(1)).^2+(he3_txy(:,3)-he3_cen(2)).^2),'sigma',sigmar,'edges',rad_bins);
            % 
            % bin_centres_he4 = count_hist_he4.bin.centers;
            % bin_centres_he3 = count_hist_he3.bin.centers;
            % dbin_r = bin_centres_he4(2)-bin_centres_he4(1);
            % 
            % count_hist_he4.count_rate.smooth = count_hist_he4.counts.smooth./(2.*pi*(bin_centres_he4*dbin_r+dbin_r^2/2));
            % count_hist_he3.count_rate.smooth = count_hist_he3.counts.smooth./(2.*pi*(bin_centres_he3*dbin_r+dbin_r^2/2));



            % flux_he4{kk} = count_hist_he4.count_rate.smooth;
            % flux_he3{kk} = count_hist_he3.count_rate.smooth;


        end
        he_dist.flux_he4{kk} = count_hist_he4.count_rate.smooth;
        he_dist.flux_he3{kk} = count_hist_he3.count_rate.smooth;
        he_dist.bin_centres_he4{kk} = count_hist_he4.bin.centers;
        he_dist.bin_centres_he3{kk} = count_hist_he3.bin.centers;
        %record atom number
        he_dist.N_he4(kk) = size(he4_txy(:,1),1);
        he_dist.N_he3(kk) = size(he3_txy(:,1),1);

        he_dist.he4_txy{kk} = he4_txy;
        he_dist.he3_txy{kk} = he3_txy;
        
        he_dist.dir_list{kk} = data.dir{ii};

        kk = kk+1;
    end
end
%%


%%
function out_folders = find_files(path,folder_list)
data_repo = dir(path);
dfolders = data_repo([data_repo(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
name_list = {dfolders.name};
types = ["atom_num","flip","dld"];
for ii = 1:length(name_list)
    if contains(name_list{ii},types)
        for jj = 1:length(types)
            if contains(name_list{ii},types(jj))
                indx = jj + 1;
            end
        end
    else
        indx = 1;
    end
    current_folder = fullfile(path,name_list{ii});
    if sum(contains({dir(current_folder).name},'.txt')&contains({dir(current_folder).name},'d'))>0
        folder_list{indx} = [folder_list{indx},current_folder];
    end
    folder_list = find_files(current_folder,folder_list);
end
out_folders = folder_list;
end

