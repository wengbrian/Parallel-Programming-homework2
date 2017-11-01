echo origin
time srun -n 2 ./ms_seq 4 -2 2 -2 2 400 400 out.png
echo mpi static
time srun -n 2 ./ms_mpi_static 4 -2 2 -2 2 400 400 out1.png
hw2-diff out.png out1.png
echo mpi dynamic
time srun -n 2 ./ms_mpi_dynamic 4 -2 2 -2 2 400 400 out2.png
hw2-diff out.png out2.png
