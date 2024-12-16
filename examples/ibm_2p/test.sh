./clean.sh
clear
./compile.sh

#sim
L=1024
t_out=0.01
t_end=020.00001
femax=1e-3
uemax=1.001
#drop
max_level=15
xd0=-800.
yd0=0.
rd0=0.5
#sphere
max_level_ibm=15
xcy=-300.
rcy=5.

#parameters
u0=37.037
rhog=1.2e-3
mul=3.70e-03
mug=6.67e-05

#run
#./impact_serial $max_level $L $u0 $t_out  $t_end $xd0 $rhog $mul $mug $femax $uemax $yd0 $xcy $rcy $max_level_ibm $rd0
#./impact_serial $max_level $L $u0 $t_out  $t_end $xd0 $rhog $mul $mug $femax $uemax $yd0 $xcy $rcy $max_level_ibm $rd0
#export HWLOC_COMPONENTS=-gl
mpiexec -np 16 ./impact_mpi $max_level $L $u0 $t_out  $t_end $xd0 $rhog $mul $mug $femax $uemax $yd0 $xcy $rcy $max_level_ibm $rd0
#mpiexec -np 16 xterm -e gdb -ex run --args ./impact_mpi $max_level $L $u0 $t_out  $t_end $xd0 $rhog $mul $mug $femax $uemax $yd0 $xcy $rcy $max_level_ibm $rd0

