import os
os.system('make')
for width in [123, 321, 1000, 1600]:
    for height in [123, 321, 1000, 1600]:
        for p in [1,2,4,8,12]:
            cmd = './ms_seq 4 -4 4 -4 4 {} {} out.png'.format(width, height)
            os.system(cmd)
            cmd = 'srun -p batch -n {} ./ms_mpi_static 4 -4 4 -4 4 {} {} out1.png'.format(p, width, height)
            print(cmd)
            os.system(cmd)
            os.system('hw2-diff out.png out1.png') 
            cmd = 'srun -p batch -n {} ./ms_mpi_dynamic 4 -4 4 -4 4 {} {} out2.png'.format(p, width, height)
            print(cmd)
            os.system(cmd)
            os.system('hw2-diff out.png out2.png') 
            cmd = 'srun -p batch -c {} -n1 ./ms_omp 4 -4 4 -4 4 {} {} out3.png'.format(p, width, height)
            print(cmd)
            os.system(cmd)
            os.system('hw2-diff out.png out3.png') 
            cmd = 'srun -p batch -c{} -n{} ./ms_omp 4 -4 4 -4 4 {} {} out4.png'.format(12//p, p, width, height)
            print(cmd)
            os.system(cmd)
            os.system('hw2-diff out.png out4.png') 
