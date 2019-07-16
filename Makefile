
PNETCDF_PATH=/opt/parallel-netcdf-1.6.1_gnu

CC      = mpic++
CFLAGS  = -O3 -I${PNETCDF_PATH}/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

all:
	${CC} ${CFLAGS} life.cpp -o gameOfLife ${LDFLAGS}

