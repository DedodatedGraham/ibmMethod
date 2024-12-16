# Serial
#qcc -g -disable-dimensions -Wall -grid=quadtree -std=c99 -O2 impact.c -o impact_serial -L$BASILISKGL -lglutils -lfb_tiny -lm
#qcc -g  -disable-dimensions -Wall -grid=quadtree -std=c99 -O2 nofimpact.c -o impact_serial -L$BASILISKGL -lglutils -lfb_tiny -lm

# openMP
#qcc -Wall -fopenmp -grid=quadtree -std=c99 -O2 impact.c -o impact_op -lm
#qcc -Wall -fopenmp -grid=quadtree -std=c99 -O2 impact.c -o impact -I$HOME -L$HOME/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

# mpi
qcc -disable-dimensions -source -D_MPI=1 -grid=quadtree impact.c
mpicc -Wall -D_XOPEN_SOURCE=700 -std=c99 -O2 _impact.c -o impact_mpi -L$BASILISKGL -lglutils -lfb_tiny -lm
#qcc -disable-dimensions -source -D_MPI=1 -grid=quadtree nofimpact.c
#mpicc -Wall -D_XOPEN_SOURCE=700 -std=c99 -O2 _nofimpact.c -o impact_mpi -L$BASILISKGL -lglutils -lfb_tiny -lm
