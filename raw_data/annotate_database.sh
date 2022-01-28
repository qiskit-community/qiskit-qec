#!/bin/bash
ADDPROPERTY=../add_property
# Do "planar" last, as it relies on "css" and "gf4linear"
PROPERTIES="standardform css decomposable degenerate triorthogonal gf4linear planar"
dmin=0
for prop in ${PROPERTIES}
do
echo COMPUTING PROPERTY ${prop} ...
# n = 1 qubits
${ADDPROPERTY} ${prop} ${dmin} database_1_0.json
mv database_1_0.json.${prop} database_1_0.json
# n = 2 qubits
${ADDPROPERTY} ${prop} ${dmin} database_2_0.json
mv database_2_0.json.${prop} database_2_0.json
${ADDPROPERTY} ${prop} ${dmin} database_2_1.json
mv database_2_1.json.${prop} database_2_1.json
# n = 3 qubits
${ADDPROPERTY} ${prop} ${dmin} database_3_0.json
mv database_3_0.json.${prop} database_3_0.json
${ADDPROPERTY} ${prop} ${dmin} database_3_1.json
mv database_3_1.json.${prop} database_3_1.json
${ADDPROPERTY} ${prop} ${dmin} database_3_2.json
mv database_3_2.json.${prop} database_3_2.json
# n = 4 qubits
${ADDPROPERTY} ${prop} ${dmin} database_4_0.json
mv database_4_0.json.${prop} database_4_0.json
${ADDPROPERTY} ${prop} ${dmin} database_4_1.json
mv database_4_1.json.${prop} database_4_1.json
${ADDPROPERTY} ${prop} ${dmin} database_4_2.json
mv database_4_2.json.${prop} database_4_2.json
${ADDPROPERTY} ${prop} ${dmin} database_4_3.json
mv database_4_3.json.${prop} database_4_3.json
# n = 5 qubits
${ADDPROPERTY} ${prop} ${dmin} database_5_0.json
mv database_5_0.json.${prop} database_5_0.json
${ADDPROPERTY} ${prop} ${dmin} database_5_1.json
mv database_5_1.json.${prop} database_5_1.json
${ADDPROPERTY} ${prop} ${dmin} database_5_2.json
mv database_5_2.json.${prop} database_5_2.json
${ADDPROPERTY} ${prop} ${dmin} database_5_3.json
mv database_5_3.json.${prop} database_5_3.json
${ADDPROPERTY} ${prop} ${dmin} database_5_4.json
mv database_5_4.json.${prop} database_5_4.json
# n = 6 qubits
${ADDPROPERTY} ${prop} ${dmin} database_6_0.json
mv database_6_0.json.${prop} database_6_0.json
${ADDPROPERTY} ${prop} ${dmin} database_6_1.json
mv database_6_1.json.${prop} database_6_1.json
${ADDPROPERTY} ${prop} ${dmin} database_6_2.json
mv database_6_2.json.${prop} database_6_2.json
${ADDPROPERTY} ${prop} ${dmin} database_6_3.json
mv database_6_3.json.${prop} database_6_3.json
${ADDPROPERTY} ${prop} ${dmin} database_6_4.json
mv database_6_4.json.${prop} database_6_4.json
${ADDPROPERTY} ${prop} ${dmin} database_6_5.json
mv database_6_5.json.${prop} database_6_5.json
# n = 7 qubits
${ADDPROPERTY} ${prop} ${dmin} database_7_0.json
mv database_7_0.json.${prop} database_7_0.json
${ADDPROPERTY} ${prop} ${dmin} database_7_1.json
mv database_7_1.json.${prop} database_7_1.json
${ADDPROPERTY} ${prop} ${dmin} database_7_2.json
mv database_7_2.json.${prop} database_7_2.json
${ADDPROPERTY} ${prop} ${dmin} database_7_3.json
mv database_7_3.json.${prop} database_7_3.json
${ADDPROPERTY} ${prop} ${dmin} database_7_4.json
mv database_7_4.json.${prop} database_7_4.json
${ADDPROPERTY} ${prop} ${dmin} database_7_5.json
mv database_7_5.json.${prop} database_7_5.json
${ADDPROPERTY} ${prop} ${dmin} database_7_6.json
mv database_7_6.json.${prop} database_7_6.json
done
