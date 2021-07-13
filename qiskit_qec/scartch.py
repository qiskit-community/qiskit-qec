# class OrderedQubits():
#     def __init__(self, qubits:[Qubit]):
#         self.qubits = sorted(qubits, cmp=self.less)
#         print(f"sorted qubits: {qubits}")
#
#
#     @property
#     def qubits(self):
#         return self.qubits
#
#     @qubits.setter
#     def set_qubits(self, qubits):
#         self.qubits = sorted(qubits, cmp=self.less)
#
#     def add_qubit(self, qb):
#         self.qubits.append(qb)
#         self.qubits = sorted(self.qubits, cmp=self.less)
#
#
#     @property
#     def centroid(self):
#         qlen = len(self.qubits)
#         x_sum = 0
#         y_sum = 0
#         for q in self.qubits:
#             x_sum += q.x()
#             y_sum += q.y()
#         centroid = QPointF(x_sum/qlen, y_sum/qlen)
#         return centroid
#
#
#     def less(self, qba:Qubit, qbb:Qubit):
#         point_a = qba.get_center()
#         point_b = qbb.get_center()
#         if (point_a.x() - self.centroid.x() >= 0) and ( point_b.x() - self.centroid.x() < 0):
#             return True
#         if (point_a.x() - self.centroid.x() < 0) and ( point_b.x() - self.centroid.x() >= 0):
#             return False
#
#         if (point_a.x() - self.centroid.x() == 0) and ( point_b.x() - self.centroid.x() == 0):
#             if (point_a.y() - self.centroid.y() >= 0) or ( point_b.y() - self.centroid.y() >= 0):
#                 return point_a.y() > point_b.y()
#             return point_b.y() > point_a.y()
#
#
#         # compute the cross product of vectors(center -> a) x(center -> b) int
#         det = (point_a.x() - self.centroid.x()) * (point_b.y() - self.centroid.y()) - (point_b.x() - self.centroid.x()) * (point_a.y() - self.centroid.y())
#         if (det < 0):
#             return True
#         if (det > 0):
#             return False
#
#         # points a and b are on the same line from the center
#
#         # check which point is closer to the center
#         d1 = (point_a.x() - self.centroid.x()) * (point_a.x() - self.centroid.x()) + (point_a.y() - self.centroid.y()) * (point_a.y() - self.centroid.y())
#         d2 = (point_b.x() - self.centroid.x()) * (point_b.x() - self.centroid.x()) + (point_b.y() - self.centroid.y()) * (point_b.y() - self.centroid.y())
#         return d1 > d2
