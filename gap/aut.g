#
# Functions for computing automorphism groups of GF(4)-additive codes
#
LoadPackage("GUAVA");

################################################
# PhiFunction: Matrix(GF(4)) -> Matrix(GF(2))
# Applies phi to each element of the input matrix, where phi is
# defined by:
# phi(0) = [0, 0, 0]
# phi(1) = [0, 1, 1]
# phi(\omega) = [1, 0, 1]
# phi(\omega^2) = [1, 1, 0]
# If a, b are elements in GF(4), then phi(a+b) = phi(a) + phi(b).
PhiFunction := function(M)
local PhiM, dims, i, j;
dims := DimensionsMat(M);
PhiM := NullMat(dims[1], 3*dims[2], GF(2));
for i in [1..dims[1]] do
  for j in [1..dims[2]] do
    if M[i,j] = Z(4)^0 then  # 1 -> 011
      PhiM[i,3*(j-1)+2] := Z(2)^0;
      PhiM[i,3*(j-1)+3] := Z(2)^0;
    elif M[i,j] = Z(4) then  # omega -> 101
      PhiM[i,3*(j-1)+1] := Z(2)^0;
      PhiM[i,3*(j-1)+3] := Z(2)^0;
    elif M[i,j] = Z(4)^2 then  # omega^2 -> 110
      PhiM[i,3*(j-1)+1] := Z(2)^0;
      PhiM[i,3*(j-1)+2] := Z(2)^0;
    fi;
  od;
od;
return PhiM;
end;

################################################
# WholeSpaceAutGroup: Integer -> Group
# Construction the automorphism group of the whole space code,
# i.e. the code containing all length-n vectors over GF(4).
# This group is S_3 \semidirect S_n. We compute the generators
# in terms of their action on the vectors in the image of Phi.
WholeSpaceAutGroup := function(n)
local i, j, plist, generators;
generators := [];
plist := [];
for j in [0..n-1] do  # Pauli permutations
  i := 3*j+1;
  Append(generators, [(i, i+1, i+2)]);
  Append(generators, [(i, i+1)]);
od;
for j in [0..3*n-1] do  # coordinate permutations
  Append(plist, [RemInt(j+3, 3*n)+1]);
od;
Append(generators, [PermList(plist)]);
if n > 1 then
  Append(generators, [(1, 4)(2, 5)(3, 6)]);
fi;
return Group(generators);
end;

################################################
# WholeSpacePermAutGroup: Integer -> Group
# Construction the permutation automorphism group of the whole
# space code, i.e. the code containing all length-n vectors
# over GF(4). This group is S_n. We compute the generators
# in terms of their action on the vectors in the image of Phi.
WholeSpacePermAutGroup := function(n)
local i, j, plist, generators;
generators := [];
plist := [];
for j in [0..3*n-1] do  # coordinate permutations
  Append(plist, [RemInt(j+3, 3*n)+1]);
od;
Append(generators, [PermList(plist)]);
if n > 1 then
  Append(generators, [(1, 4)(2, 5)(3, 6)]);
fi;
return Group(generators);
end;

################################################
# WholeSpaceS3nAutGroup: Integer -> Group
# Construction the subgroup of the automorphism group of the
# whole space code corresponding to relabeling the non-zero
# elements of GF(4). This group is S_3^n. We compute the
# generators in terms of their action on the vectors in the
# image of Phi.
WholeSpaceS3nAutGroup := function(n)
local i, j, generators;
generators := [];
for j in [0..n-1] do  # Pauli permutations
  i := 3*j+1;
  Append(generators, [(i, i+1, i+2)]);
  Append(generators, [(i, i+1)]);
od;
return Group(generators);
end;

################################################
# GF4AdditiveAutGroup: Matrix(GF(4)) -> Group
# The input matrix G is an n-k by n matrix over GF(4)
# that generates the GF(4)-additive code, i.e. the code
# is all GF(2) linear combinations of the rows of G.
# Returns the automorphism group of the code.
GF4AdditiveAutGroup := function(G)
local PhiC, dims, AC, AW;
dims := DimensionsMat(G);
PhiC := GeneratorMatCode(PhiFunction(G), GF(2));  # a linear code
AC := AutomorphismGroup(PhiC);
AW := WholeSpaceAutGroup(dims[2]);
return Intersection(AC, AW);
end;

################################################
# GF4AdditivePermAutGroup: Matrix(GF(4)) -> Group
# The input matrix G is an n-k by n matrix over GF(4)
# that generates the GF(4)-additive code, i.e. the code
# is all GF(2) linear combinations of the rows of G.
# Returns the permutation automorphism group of the code.
GF4AdditivePermAutGroup := function(G)
local PhiC, dims, AC, AW;
dims := DimensionsMat(G);
PhiC := GeneratorMatCode(PhiFunction(G), GF(2));  # a linear code
AC := AutomorphismGroup(PhiC);
AW := WholeSpacePermAutGroup(dims[2]);
return Intersection(AC, AW);
end;

################################################
# GF4AdditiveLocalAutGroup: Matrix(GF(4)) -> Group
# The input matrix G is an n-k by n matrix over GF(4)
# that generates the GF(4)-additive code, i.e. the code
# is all GF(2) linear combinations of the rows of G.
# Returns the subgroup of the automorphism group of the
# code that relabels non-zero elements of GF(4) on
# each coordinate. The elements of this group correspond
# to bitwise Clifford operations that perserve the code
# space of a stabilizer code.
GF4AdditiveLocalAutGroup := function(G)
local PhiC, dims, AC, AW;
dims := DimensionsMat(G);
PhiC := GeneratorMatCode(PhiFunction(G), GF(2));  # a linear code
AC := AutomorphismGroup(PhiC);
AW := WholeSpaceS3nAutGroup(dims[2]);
return Intersection(AC, AW);
end;

################################################
# StringsToGF4: List(string) -> Matrix(GF(4))
# Convert a list of strings of 'x', 'y', and 'z'
# characters representing the Pauli operators into
# a matrix over GF(4) whose rows correspond to
# the strings. Assume that Length(stringlist[i])
# is the same for all i.
StringsToGF4 := function(stringlist)
local M, i, j;
M := NullMat(Length(stringlist), Length(stringlist[1]), GF(4));
for i in [1..Length(stringlist)] do
  for j in [1..Length(stringlist[1])] do
    if stringlist[i][j] = 'x' then
      M[i][j] := Z(4)^0;
    elif stringlist[i][j] = 'z' then
      M[i][j] := Z(4);
    elif stringlist[i][j] = 'y' then
      M[i][j] := Z(4)^2;
    fi;
  od;
od;
return M;
end;

################################################
# reformat: (Group, Integer) -> Record
# Produce a more human-readable representation of the
# generators of (a subgroup of) the automorphism group
# of a GF(4)-additive code. The record contains a list
# of coordinate permutations and a corresponding list of
# local relabelings of the non-zero elements of GF(4).
# These can be interpreted as pairs of coordinate permutations
# preceeded by single-qubit Clifford gates that, taken together,
# map the stabilizer code space to itself. In the returned
# record, the Clifford gates are encoded as:
# 'i': identity
# 'r': 1 -> \omega -> \omega^2 (x -> z -> y)
# 'R': 1 -> \omega^2 -> \omega (x -> y -> z)
# 'V': \omega <-> \omega^2 (y <-> z)
# 'S': 1 <-> \omega^2 (x <-> y)
# 'H': 1 <-> \omega (x <-> z)
reformat := function(A, n)
  local autgens, g, img, perm, sing, sing2, tup, permlist, lclist, i;
  autgens := GeneratorsOfGroup(A);
  permlist := [];
  lclist := [];
  for g in autgens do
    img := OnTuples([1..3*n], g);
    perm := [1..n];
    sing := [];
    for i in [1..n] do
      perm[QuoInt(img[3*(i-1)+1]-1, 3)+1] := i;
      tup := [RemInt(img[3*(i-1)+1]-1, 3),
              RemInt(img[3*(i-1)+2]-1, 3),
              RemInt(img[3*(i-1)+3]-1, 3)];
      if tup = [0, 1, 2] then
        Append(sing, "I");  # *1
      elif tup = [1, 2, 0] then
        Append(sing, "R");  # *a^2 = *Z(4)^2 = b
      elif tup = [2, 0, 1] then
        Append(sing, "r");  # *a = *Z(4)
      elif tup = [1, 0, 2] then
        Append(sing, "H");  # exchange 1 and a, fix a^2
      elif tup = [0, 2, 1] then
        Append(sing, "V");  # exchange a and a^2, fix 1
      elif tup = [2, 1, 0] then
        Append(sing, "S");  # exchange 1 and a^2, fix a
      fi;
    od;
    # Apply perm^-1 to sing
    sing2 := [];
    for i in [1..n] do
      Append(sing2, [sing[perm[i]]]);
    od;
    Add(permlist, perm);
    Add(lclist, sing2);
  od;
  return rec(permutations:=permlist, locals:=lclist);
end;
