
#PNETCDF_PATH=/opt/parallel-netcdf-1.6.1_gnu
PNETCDF_PATH=/opt/cray/pe/parallel-netcdf/1.8.1.3/INTEL/16.0

CC      = CC
CFLAGS  = -O3 -profile-loops=all -I${PNETCDF_PATH}/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf -upthread_attr_getstacksize -upthread_attr_setstack -upthread_getattr_np -upthread_attr_getstack

all:
	${CC} ${CFLAGS} life.cpp -o gameOfLife ${LDFLAGS}

