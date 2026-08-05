#!/bin/bash
#SBATCH --job-name=Find_motifs_thresh_all-MW
#SBATCH --output=find_motifs_thresh_%a.out
#SBATCH -t 12:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12

#SBATCH --array=1-5

ml load matlab
matlab -nodisplay -nosplash -nodesktop -r "\
load('FlexMF_params_all.mat'); \
disp(['Task ID: ', num2str($SLURM_ARRAY_TASK_ID)]);\
input_folder_names = {'20240618MH041724', '20240602FH040624', '20240326ML102623', '20250127ML111624', '20250114FL111424'};\
input_folder_name = input_folder_names{$SLURM_ARRAY_TASK_ID};\
disp(['Input Folder Name: ', input_folder_name]);\
p = parameters{input_folder_name, :};\
disp(['Parameters: ', num2str(p)]);\
if isempty(regexp(input_folder_name, '.+ML102623', 'once')),\
    L = 100;\
else,\
    L = 50;\
end;\
disp(['L: ', num2str(L)]);\
find_motifs(input_folder_name, 'lambda', p(1), 'alpha_W', p(2), 'alpha_H', p(3), 'L', L, 'calcium_thresh', 1.5, 'local', false);\
find_motifs(input_folder_name, 'lambda', p(1), 'alpha_W', p(2), 'alpha_H', p(3), 'L', L, 'calcium_thresh', 2.5, 'local', false);\
exit;\
"
