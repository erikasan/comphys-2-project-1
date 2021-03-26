
# Example runs
This folder contains the files that are madeafter running the following commands:
 ```bash
 mkdir example_runs
./vmc 3 10 10000 0.5 0.1 100000 100 HO IMP no ./example_runs/
./vmc 3 10 10000 0.85 0.1 100000 100 HO VMC no ./example_runs/
./vmc 3 10 262144 0.4975 0.1 100000 100 EO IMP energyanddistribution ./example_runs/
./vmc_parallel 3 10 262144 0.4975 0.1 100000 100 EO IMP energiesparallelrun 4 ./example_runs/
./gradientdescent 3 10 10000 0.45 0.1 100000 100 HO IMP gradientdescent ./example_runs/
./gradientdescent_parallel 3 10 10000 0.55 0.1 100000 100 HO IMP gradientdescent 4 ./example_runs/
 ```
 
 Observe that the file "alphas.csv" is directly generated from the main file, while the other files are created by the System. The filename "no" creates no files.
