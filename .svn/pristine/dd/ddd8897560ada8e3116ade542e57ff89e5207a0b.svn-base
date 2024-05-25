FC = gfortran
SRC = 	fortran/common.f90 \
		fortran/kdbinsert.f90
FMOD = 	fortran/fhelper.f90
OPT = -O2

main : fortran_modules
	f2py -c --opt="${OPT}" -m fmodules ${SRC} ${FMOD}

fortran_modules : ${FMOD}
	${FC} ${OPT} -c ${FMOD}

clean :
	rm -f fmodules.*so
	rm -f *.o
	rm -f *.mod
