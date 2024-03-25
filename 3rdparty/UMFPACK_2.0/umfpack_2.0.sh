#!/bin/bash
#
mkdir src
cd src
../f77split ../umfpack_2.0.f
#
for FILE in `ls -1 *.f`;
do
  gfortran -c $FILE -g -fdefault-real-8 -O2
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
#
ar qc libumfpack_2.0.a *.o
rm *.o
#
mv libumfpack_2.0.a ../
cd ..
#
echo "Library installed"
