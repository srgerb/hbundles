sbatch -a 1-$(cat array_tasks.list | wc -l) sbatch_array_job.sh
