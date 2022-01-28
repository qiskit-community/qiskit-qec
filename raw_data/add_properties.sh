#!/bin/bash
ADDPROPERTY=../add_property
# Do "planar" last, as it relies on "css" and "gf4linear"
PROPERTIES="standardform decomposable degenerate planar"
#PROPERTIES="standardform decomposable degenerate css triorthogonal gf4linear planar"
dmin=0
for prop in ${PROPERTIES}
do
echo COMPUTING PROPERTY ${prop} ...
${ADDPROPERTY} ${prop} ${dmin} $1 
mv $1.${prop} $1
sleep 5
done
