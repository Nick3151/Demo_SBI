#!/bin/bash
#SBATCH --job-name=Analyze_motifs_shuffled-MW
#SBATCH --output=analyze_motifs_shuffled_across_trials_%a.out
#SBATCH -t 24:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12

#SBATCH --array=1-40

ml load matlab
matlab -nodisplay -nosplash -nodesktop -r "\
load('FlexMF_params_all.mat'); \
disp(['Task ID: ', num2str($SLURM_ARRAY_TASK_ID)]);\
input_folder_names = {'20240618MH041724', '20240624MH041724', '20240702MH041724', '20240530FH040624', '20240602FH040624', '20240709MH041724','20240611FH040624', '20240307ML102623', '20240310ML102623', '20240326ML102623','20240328ML102623', '20240405ML102623', '20241011ML072024', '20241016ML072024','20241101ML072024', '20241111ML072024', '20241123ML072024', '20241205ML072024', '20241011FL073024', '20250101MH111624', '20250102MH111624', '20250104MH111624', '20250105MH111624', '20250108MH111624', '20250110MH111624', '20250109ML111624', '20250113ML111624', '20250116ML111624', '20250119ML111624', '20250120ML111624', '20250122ML111624', '20250114FL111424', '20250116FL111424', '20250118FL111424', '20250120FL111424', '20250117ML111624', '20250123ML111624', '20250126ML111624', '20250127ML111624', '20250130ML111624'};\
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
save_file = 'motifs_shuffled_across_trials.mat';\
find_motifs_shuffle(input_folder_name, 'across_trials', save_file, 'lambda', p(1), 'alpha_W', p(2), 'alpha_H', p(3), 'L', L, 'local', false, 'seed', 1);\

exit;\
"
