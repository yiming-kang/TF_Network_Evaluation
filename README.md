# TF Network Evaluation

This repo analysis the performance of network inference output against benchmark datasets: ChIP-network and PWM-network. It outputs evaluation figures as rank-precision plotting format. It contains two steps:

#### 1. Analyze performance
```
sbatch --mail-type=FAIL,END --mail-user=<email> evaluate_network_w_nonbinary_benchmark.sh <network_file> <regulators_file> <genes_file> 32 20 <output_directory>
```

#### 2. Make evaluation plots
```
python plot_evaluations.py --network_evals <network_eval_list> --network_labels <label_list> --random_eval_dir <random_eval_dir> --figure_file_suffix <figure_suffix> --step 1600 --num_regulators 320
```

