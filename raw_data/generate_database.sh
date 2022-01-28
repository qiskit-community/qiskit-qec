#!/bin/bash
CLASSIFY=../classify
CLASSIFY2=../classify2
GRAPHS=../third_party/graphs
# n = 1 qubits
${CLASSIFY} database_1_0.json 0 ${GRAPHS}/n1/n1.dimacs
# n = 2 qubits
${CLASSIFY} database_2_0.json 0 ${GRAPHS}/n2/n2-d2.dimacs ${GRAPHS}/n2/decom/n2-d2-11.dimacs
${CLASSIFY} database_2_1.json 1 ${GRAPHS}/n2/n2-d2.dimacs ${GRAPHS}/n2/decom/n2-d2-11.dimacs
# n = 3 qubits
${CLASSIFY} database_3_0.json 0 ${GRAPHS}/n3/n3-d2.dimacs ${GRAPHS}/n3/decom/n3-1*.dimacs
${CLASSIFY} database_3_1.json 1 ${GRAPHS}/n3/n3-d2.dimacs ${GRAPHS}/n3/decom/n3-1*.dimacs
${CLASSIFY} database_3_2.json 2 ${GRAPHS}/n3/n3-d2.dimacs ${GRAPHS}/n3/decom/n3-1*.dimacs
# n = 4 qubits
${CLASSIFY} database_4_0.json 0 ${GRAPHS}/n4/*.dimacs ${GRAPHS}/n4/decom/n4-*.dimacs
${CLASSIFY} database_4_1.json 1 ${GRAPHS}/n4/*.dimacs ${GRAPHS}/n4/decom/n4-*.dimacs
${CLASSIFY} database_4_2.json 2 ${GRAPHS}/n4/*.dimacs ${GRAPHS}/n4/decom/n4-*.dimacs
${CLASSIFY} database_4_3.json 3 ${GRAPHS}/n4/*.dimacs ${GRAPHS}/n4/decom/n4-*.dimacs
# n = 5 qubits
${CLASSIFY} database_5_0.json 0 ${GRAPHS}/n5/n5-*.dimacs ${GRAPHS}/n5/decom/n5-*.dimacs
${CLASSIFY} database_5_1.json 1 ${GRAPHS}/n5/n5-*.dimacs ${GRAPHS}/n5/decom/n5-*.dimacs
${CLASSIFY} database_5_2.json 2 ${GRAPHS}/n5/n5-*.dimacs ${GRAPHS}/n5/decom/n5-*.dimacs
${CLASSIFY} database_5_3.json 3 ${GRAPHS}/n5/n5-*.dimacs ${GRAPHS}/n5/decom/n5-*.dimacs
${CLASSIFY} database_5_4.json 4 ${GRAPHS}/n5/n5-*.dimacs ${GRAPHS}/n5/decom/n5-*.dimacs
# n = 6 qubits
${CLASSIFY} database_6_0.json 0 ${GRAPHS}/n6/n6-*.dimacs ${GRAPHS}/n6/decom/n6-*.dimacs
${CLASSIFY} database_6_1.json 1 ${GRAPHS}/n6/n6-*.dimacs ${GRAPHS}/n6/decom/n6-*.dimacs
${CLASSIFY} database_6_2.json 2 ${GRAPHS}/n6/n6-*.dimacs ${GRAPHS}/n6/decom/n6-*.dimacs
${CLASSIFY} database_6_3.json 3 ${GRAPHS}/n6/n6-*.dimacs ${GRAPHS}/n6/decom/n6-*.dimacs
${CLASSIFY} database_6_4.json 4 ${GRAPHS}/n6/n6-*.dimacs ${GRAPHS}/n6/decom/n6-*.dimacs
${CLASSIFY2} database_6_5.json 6 5
# n = 7 qubits
${CLASSIFY} database_7_0.json 0 ${GRAPHS}/n7/n7-*.dimacs ${GRAPHS}/n7/decom/n7-*.dimacs
${CLASSIFY} database_7_1.json 1 ${GRAPHS}/n7/n7-*.dimacs ${GRAPHS}/n7/decom/n7-*.dimacs
${CLASSIFY} database_7_2.json 2 ${GRAPHS}/n7/n7-*.dimacs ${GRAPHS}/n7/decom/n7-*.dimacs
${CLASSIFY} database_7_3.json 3 ${GRAPHS}/n7/n7-*.dimacs ${GRAPHS}/n7/decom/n7-*.dimacs
${CLASSIFY} database_7_4.json 4 ${GRAPHS}/n7/n7-*.dimacs ${GRAPHS}/n7/decom/n7-*.dimacs
${CLASSIFY2} database_7_5.json 7 5
${CLASSIFY2} database_7_6.json 7 6
# n = 8 qubits
${CLASSIFY} database_8_0.json 0 ${GRAPHS}/n8/n8-*.dimacs ${GRAPHS}/n8/decom/n8-*.dimacs
${CLASSIFY} database_8_1.json 1 ${GRAPHS}/n8/n8-*.dimacs ${GRAPHS}/n8/decom/n8-*.dimacs
${CLASSIFY} database_8_2.json 2 ${GRAPHS}/n8/n8-*.dimacs ${GRAPHS}/n8/decom/n8-*.dimacs
${CLASSIFY} database_8_3.json 3 ${GRAPHS}/n8/n8-*.dimacs ${GRAPHS}/n8/decom/n8-*.dimacs
# missing _8_4
# missing _8_5
${CLASSIFY2} database_8_6.json 8 6
${CLASSIFY2} database_8_7.json 8 7

