#Example runs
This folder contains the files that are madeafter running the following commands:
 ```bash
./vmc 3 10 10000 0.5 0.1 100000 100 HO IMP no
./vmc 3 10 10000 0.85 0.1 100000 100 HO VMC no
./vmc 3 10 262144 0.4975 0.1 100000 100 EO IMP energyanddistribution
./vmc_parallel 3 10 262144 0.4975 0.1 100000 100 EO IMP energiesparallelrun 4
./gradientdescent 3 10 10000 0.45 0.1 100000 100 HO IMP gradientdescent
./gradientdescent_parallel 3 10 10000 0.55 0.1 100000 100 HO IMP gradientdescent 4
 ```
 
 Observe that the file "alphas.csv" is directly generated from the main file, while the other files are created by the System. The filename "no" creates no files.
