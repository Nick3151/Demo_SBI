#!/bin/bash
#SBATCH --job-name=Analyze_motifs_shuffled_within-MW
#SBATCH --output=analyze_motifs_shuffled_within_trials_%a.out
#SBATCH -t 24:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12

#SBATCH --array=1-100

ml load matlab
matlab -nodisplay -nosplash -nodesktop -r "\
load('FlexMF_params_all.mat'); \
task_id = $SLURM_ARRAY_TASK_ID;\
disp(['Task ID: ', num2str(task_id)]);\
input_folder_names = {'20240618MH041724', '20240602FH040624', '20240326ML102623', '20250127ML111624', '20250114FL111424'};\
nSessions = length(input_folder_names);\
n = mod(task_id, nSessions);\
if n==0,\
    n = 5;\
end;\
input_folder_name = input_folder_names{n};\
disp(['Input Folder Name: ', input_folder_name]);\
p = parameters{input_folder_name, :};\
disp(['Parameters: ', num2str(p)]);\
if isempty(regexp(input_folder_name, '.+ML102623', 'once')),\
    L = 100;\
else,\
    L = 50;\
end;\
disp(['L: ', num2str(L)]);\
seed = ceil(task_id/nSessions);\
save_file = sprintf('motifs_shuffled_within_trials_%02d.mat', seed);\
find_motifs_shuffle(input_folder_name, 'within_trials', save_file, 'lambda', p(1), 'alpha_W', p(2), 'alpha_H', p(3), 'L', L, 'local', false, 'seed', seed);\
exit;\
"
