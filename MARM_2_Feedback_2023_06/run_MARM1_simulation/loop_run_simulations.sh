for i in `seq 0 0`; do
  sbatch single_run_simulation.sh parameter_sets.csv $i
done

for i in `seq 1 99`; do
  sbatch single_run_reduced_simulation.sh parameter_sets.csv $i
done
