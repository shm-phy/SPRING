
FC=gfortran

#sed -i "1s/^/FC=$FC\n/" Makefile
COARR_LIB=/opt/opencoarrays/lib

echo $(subst /,\/,${COARR_LIB})


