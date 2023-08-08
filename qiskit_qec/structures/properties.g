LoadPackage("Gauss");


# Create a GF(2) Lambda Matrix
# [0 I_(n/2)]
# [I_(n/2) 0]
# where I_n is a nxn identity matrix
CreateLambda := function(n)
  local i, lambda;

  lambda := NullMat(n,n);

  for i in [1..n/2] do
    lambda[i][n/2+i] := 1;
  od;
  for i in [n/2+1..n] do
    lambda[i][i-n/2] := 1;
  od;

  return lambda*One(GF(2));
end;

# Given an set of generators S=<g_1,g_2,...,g_k> find a hyperbolic element g in P_n to g_i
# So g commutes with every generator g_j for i different from i but anticommutes with g_i
# Reference Prop 10.4. Nielsen and Chuang
# Assumption is that S contains order two elements all of which commute
# Input S as a check matrix kx2n
FindHyperbolicPartner := function(S,index)
  local i,lambda, slambda, rref, ek, rek, pivot, ans;

  lambda := CreateLambda(Size(S[1]));

  slambda := S*lambda;

  rref := EchelonMatTransformation(slambda);

  ek := 0*[1..Size(S)];
  ek[index] := 1;
  ek := ek*One(GF(2));

  rek := rref.coeffs*ek;

  pivot := 1;
  ans := 0*[1..Size(S[1])];
  for i in [1..Size(S[1])] do
    if rref.heads[i] > 0 then
      ans[i] := rek[pivot];
      pivot := pivot + 1;
    fi;
  od;

  return ans*One(GF(2));
end;


MakeStab := function(list, flag)
  local u, out;

  out := [];

  for u in list do
    if flag = "X" then
      Add(out, [u,[]]);
    else
      Add(out, [[],u]);
    fi;
  od;
  return out;
end;

AddElementsOf := function(listo, listi)
  local u;

  for u in listi do 
    Add(listo, u);
  od;
end;

SqStab := function(a,b)
  return [a,a+1,b,b+1];
end;


# Create sympletic matrix from compress format
CreateGenMatrix := function(stab,nqubits) 
  local s, u,v,k,mat,row;
  mat:=[];
  for s in stab do  
    u := s[1];
    v := s[2];
    row := 0*[1..2*nqubits];
    for k in u do
      row[k] := 1;
    od;
    for k in v do
      row[k+nqubits] := 1;
    od;
    Add(mat,row);
  od;
  return mat;
end;


# Given a matrix M agument it with a identity matrix
# [ M ]
# [ I ]
# flag = true to augment GF(2) entries not 0s and 1s
AugmentMat := function(mat,flag) 
  local u, row, nmat, k;

  nmat := StructuralCopy(mat);

  for k in [1..Size(nmat[1])] do
    row := 0*[1..Size(nmat[1])];
    row[k] := 1;
    if flag then
      row := row * One(GF(2));
    fi;
    Add(nmat,row);
  od;
  return nmat;
end;

# Compute the Sympletic product (symplectic bilinear form) or two 
# GF(2) vectors and return integer equivalent values 
SymplecticProduct := function(avec, bvec)
  local p1,p2,k,n;
  if  not Size(avec) = Size(bvec) then
    Error("Vectors are not of the same dimension");
  elif not IsEvenInt(Size(avec)) then 
    Error("Vectors are not of even dimension");
  fi;

  n :=  Size(avec)/2;

  p1:=avec[1]*bvec[n+1];
  p2:=avec[n+1]*bvec[1];

  for k in [2..n] do
    p1 := p1 + avec[k]*bvec[n+k];
    p2 := p2 + avec[n+k]*bvec[k];
  od;

  if IsOne(p1+p2) then 
    return 1;
  else 
    return 0;
  fi;
end;

# Make the operator elem commute with the hyperbolic pair Xe, Ze
MakeElementCommute := function(elem, Xe, Ze)
  if SymplecticProduct(elem, Xe) = 1 then
    elem := elem + Ze;
    # elem := GroupMultiplication(elem, Ze);
  fi;
  if SymplecticProduct(elem, Ze) = 1 then
    elem := elem + Xe;
    # elem := GroupMultiplication(elem, Xe);
  fi;
  return elem;
end;

# Make the operator elem commute with the hyperbolic pairs Xe(xrange) and Ze(zrange)
MakeElementCommuteGeneral := function(elem, Xe, xrange, Ze, zrange)
  local  i, elem_mod;

  elem_mod := StructuralCopy(elem);

  for i in [1..Size(xrange)] do
    elem_mod := MakeElementCommute(elem_mod, Xe[xrange[i]], Ze[zrange[i]]);
  od;

  return elem_mod;
end;

# Make the elements in S(start..finish) commute with the hyperbolic pair Xe, Ze
MakeElementsCommute := function(S, start, finish, Xe, Ze)
  local k;

  for k in [start..finish] do
    S[k] := MakeElementCommute(S[k], Xe, Ze);
  od;

  return S;
end;

# Make the elements in S(range) commute with the hyperbolic pair Xe, Ze
MakeElementsCommuteGeneral := function(S, range, Xe, Ze)
  local k;

  for k in range do
    S[k] := MakeElementCommute(S[k], Xe, Ze);
  od;

  return S;
end;

# Find generator in S that does not commute with operator el
FindNonCommunativePartner := function(S, el)
  local i;

  for i in [1..Size(S)] do
    if SymplecticProduct(S[i],el)=1 then
      return [S[i],i];
    fi;
  od;
  return [];
end;


# Apply the modified GramSchmidt process to a input matrix G_in
# Optional arg enables one to start at a specific position in the process
GramSchmidt := function(G_in, arg...)
  local G, Xe, Ze, center, size, elem, res, start, m, k;
  G:=StructuralCopy(G_in);

  if Length(arg)=1 then
    start := arg[1];
  elif Length(arg)=0 then
    start := 1;
  fi;

  m := (start-1)/2;

  Xe:=[];
  Ze:=[];
  center := [];

  for k in [1..m] do
    Xe[k] := G[k];
    Ze[k] := G[m+k];
  od;

  for k in [1..2*m] do
    Remove(G, 1);
  od;

  size := Size(G);

  while size > 0 do
    elem := G[1];
    Remove(G,1);
    size := size - 1;
    res := FindNonCommunativePartner(G, elem);
    if Size(res) = 0 then
      Add(center, elem);
    else
      Add(Xe, elem);
      Add(Ze, res[1]);
      Remove(G, res[2]);
      size := size - 1;
      G := MakeElementsCommute(G, 1, Size(G), elem, res[1]);
    fi;
  od;
  return rec(center := center, hyperX :=  Xe, hyperZ := Ze);
end;


GetCenter := function(G)
  return GramSchmidt(G).center;
end;

# Extract rows from D based on extact_list
ExtractElements := function(D, extract_list)
  local k, G;

  G := [];
  for k in [1..Size(extract_list)] do
    if extract_list[k] > 0 then
      Add(G, D[k]);
    fi;
  od;

  return G;
end;

# Extend the basis given by G to lin. ind. set for P_n
ExtendBasis := function(G)
  local G_aug, G_augt, G_augt_rref, G_ext;

  G_aug := AugmentMat(G, true);
  G_augt := TransposedMat(G_aug);
  G_augt_rref := EchelonMat(G_augt);

  G_ext := ExtractElements(G_aug, G_augt_rref.heads);

  return G_ext;
end;

# Delete elements from a hyperbolic basis (move to other elements to center as appropriate)
RemoveHyperbolicBasisElements := function(G, component, range)
  local i, lrange, G_mod, source, sink, remove, item;

  lrange := StructuralCopy(range);
  G_mod := StructuralCopy(G);

  if component = "HyperX" then
    source := G_mod.hyperZ;
    remove := [G_mod.hyperX, G_mod.hyperZ];
    sink := G_mod.center;
  elif component = "HyperZ" then
    source := G_mod.hyperX;
    remove := [G_mod.hyperX, G_mod.hyperZ];
    sink := G_mod.center;
  elif component = "Center" then
    source := [];
    remove := [G_mod.center];
    sink := [];
  fi;

  Sort(List(lrange));

  if Size(source) > 0 then
    for i in lrange do
      Add(sink, source[i]);
    od;
  fi;

  lrange := Reversed(lrange);

  for item in remove do
    for i in lrange do
      Remove(item, i);
    od;
  od;

  return G_mod;
end;


# Join (stack) two matrices together
# A, B, C, ... -> [A]
#                 [B]
#                 [.]
#                 [.]
#                 [.]
CombineComponents := function(args...)
  local k, H, object;

  H:=[];

  for object in args do
    for k in [1..Size(object)] do
      Add(H, object[k]);
    od;
  od;

  return H;
end;

# Make element v commute with hyperbolic pairs A(range) B(range)
MakeElementCommuteRestricted := function(v, A, B, range)
  local i, v_mod;

  v_mod := StructuralCopy(v);

  for i in range do
    MakeElementCommute(v_mod, A[i], B[i]);
  od;

  return v_mod;
end;



# A will be modified to commute across B (minus hyperbolic partner)
MakeHyperbolicPairsCommute := function(A,B)
  local A_mod, i;

  A_mod := StructuralCopy(A);

  for i in [1..Size(A_mod)] do
    A_mod[i] := MakeElementCommuteRestricted(A_mod[i], A_mod, B, [(i+1)..Size(B)]);
  od;

  return A_mod;
end;

MakeSetsCommuteGeneral := function(operators, orange, hyperX, xrange, hyperZ, zrange)
  local i, operators_mod;

  operators_mod := StructuralCopy(operators);

  for i in [1..Size(xrange)] do
    MakeElementsCommuteGeneral(operators_mod, orange, hyperX[xrange[i]], hyperZ[zrange[i]]);
  od;

  return operators_mod;
end;



MakeSetsCommute := function(operators, hyperX, hyperZ)
  local i, operators_mod;

  operators_mod := StructuralCopy(operators);

  for i in [1..Size(hyperX)] do
    MakeElementsCommute(operators_mod,1, Size(operators_mod), hyperX[i], hyperZ[i]);
  od;

  return operators_mod;
end;

ExtendHyperbolicBasis := function(H)
  local H_mod, i, hop;

  H_mod := StructuralCopy(H);

  while Size(H_mod.center)>0 do
    hop := FindHyperbolicPartner(H_mod.center, 1);
    hop := MakeElementCommuteGeneral(hop, H_mod.hyperX, [1..Size(H_mod.hyperX)], H_mod.hyperZ, [1..Size(H_mod.hyperZ)]);
    Add(H_mod.hyperX, H_mod.center[1]);
    Add(H_mod.hyperZ, hop);
    Remove(H_mod.center, 1);
  od;

  return H_mod;
end;



# ----------------------------------------

IsotropicHyperboicForm := function(G)
  return GramSchmidt(G);
end;

IsotropicHyperboicBasis := function(G, arg...)
  local G_int, center_size, gauge_size, G_enb, G_ext, orange, xrange, zrange, G_ext_com, G_hyp, range,P_n, G_mod;

  if Size(arg)>0 and arg[1]=false then
    G_mod := GramSchmidt(G);
  else
    G_mod := G;
  fi;
    
  center_size := Size(G_mod.center);
  gauge_size := Size(G_mod.hyperX);

  G_enb := ExtendHyperbolicBasis(G);
  G_ext := ExtendBasis(CombineComponents(G_enb.hyperX, G_enb.hyperZ));

  orange := [2*Size(G_enb.hyperX)+1..Size(G_ext)];
  xrange := [1..Size(G_enb.hyperX)];
  zrange := [Size(G_enb.hyperX)+1..2*Size(G_enb.hyperX)];

  G_ext_com := MakeSetsCommuteGeneral(G_ext, orange, G_ext, xrange, G_ext, zrange);
  P_n := GramSchmidt(G_ext_com, 2*Size(G_enb.hyperX)+1);

  return P_n;
end;

# Note order of non center generators may change
GaugeGroupCenter := function(G, arg...)
  local i,j,H, flag, irange, break_flag;

  H:=[];

  if Size(arg)=0 or arg[1]=false then
    return GramSchmidt(G).center;
  else
    irange := Reversed([1..Size(G)]);
    for i in irange do
      break_flag := 0;
      for j in [1..Size(G)] do
        if SymplecticProduct(G[i], G[j])=1 then
          break_flag := 1;
          break;
        fi;
      od;
      if break_flag = 0 then
        Add(H,G[i],1); # Add G[i] at the front of the list H (G[i] in center)
      else
        Add(H,G[i]); # Add G[i] to the end of the list H (G[i] not in center)
      fi;
    od;
    return GramSchmidt(H).center;
  fi;
end;


# G, qec_info, true/false 
# true -> maintain form
# false -> ignore form
# default = false

GaugeGroupNormalizer := function(G, qec_info, arg...)
local G_int, G_ext, G_hyp, G_ext_hyp, G_ext_hyb, G_tmp, distance_of_center, N_G, gauge_degree,
gauge_degrees, number_of_logical_qubits,Y,range, G_ext_com, orange, xrange, zrange, G_enb, center_size;

  if qec_info.isstabilizergroup = true then
    G_tmp := G;
    distance_of_center := Rank(G);
    gauge_degree := 0;

    G_ext := ExtendBasis(G_tmp);
    G_ext_hyb := GramSchmidt(G_ext);

    N_G := RemoveHyperbolicBasisElements(G_ext_hyb, "HyperZ", [1..Size(G)]);

    number_of_logical_qubits := Size(N_G.hyperX)-gauge_degree;

    return rec(num_qubits:=number_of_logical_qubits, center_distance:=distance_of_center, gauge_degree:=gauge_degree, normalizer:=N_G);
  elif Size(arg)>0 then
    if arg[1] = true then # Maintain form
      center_size := Size(G.center);
      gauge_degree := Size(G.hyperX);
      G_enb := ExtendHyperbolicBasis(G);

      G_ext := ExtendBasis(CombineComponents(G_enb.hyperX, G_enb.hyperZ));

      orange := [2*Size(G_enb.hyperX)+1..Size(G_ext)];
      xrange := [1..Size(G_enb.hyperX)];
      zrange := [Size(G_enb.hyperX)+1..2*Size(G_enb.hyperX)];

      G_ext_com := MakeSetsCommuteGeneral(G_ext, orange, G_ext, xrange, G_ext, zrange);
      G_hyp := GramSchmidt(G_ext_com, 2*Size(G_enb.hyperX)+1);

      range := [1+gauge_degree..center_size+gauge_degree];
      N_G := RemoveHyperbolicBasisElements(G_hyp, "HyperZ", range);

      number_of_logical_qubits := Size(N_G.hyperX)-gauge_degree;

      return rec(num_qubits:=number_of_logical_qubits, center_distance:=center_size, gauge_degree:=gauge_degree, normalizer:=N_G);
    fi;
  else # Ignore form
    G_tmp := GramSchmidt(G);

    # Note: technically distance of Z(G) is distance_of_center + 1 (as iI is a needed generator)
    distance_of_center := Size(G_tmp.center);
    gauge_degree := Size(G_tmp.hyperX);

    G_ext := ExtendBasis(CombineComponents(G_tmp.hyperX, G_tmp.hyperZ, G_tmp.center));

    orange := [(2*Size(G_tmp.hyperX)+1)..Size(G_ext)];
    xrange := [1..Size(G_tmp.hyperX)];
    zrange := [Size(G_tmp.hyperX)+1..2*Size(G_tmp.hyperX)];

    G_ext_com := MakeSetsCommuteGeneral(G_ext, orange , G_ext, xrange, G_ext, zrange);
    G_ext_hyp := GramSchmidt(G_ext_com, 2*Size(G_tmp.hyperX)+1);

    range := [1+gauge_degree..distance_of_center+gauge_degree];
    N_G := RemoveHyperbolicBasisElements(G_ext_hyp, "HyperZ", range);

    number_of_logical_qubits := Size(N_G.hyperX)-gauge_degree;

    return rec(num_qubits:=number_of_logical_qubits, center_distance:=distance_of_center, gauge_degree:=gauge_degree, normalizer:=N_G);
  fi;
end;

# Some Data examples (using different creation methods)

# PRP X 7,021029 (2017) Example -5 with corners corrected to give a single boundary

stabZPRP5 := [
  [1,2], [3,4], [5,6],
  [1,7,8], [2,3,9,10], [4,5,11,12], [6,13],
  [7,14], [8,9,15,16], [10,11,17,18], [12,13,19,20],
  [14,15,21,22], [16,17,23], [18,19,24,25], [20,26],
  [21,27], [22,23,28,29], [24,30,31], [25,26,32,33],
  [27,28,34,35], [29,30,36,37], [31,32,38,39], [33,40],
  [34,41], [35,36,42,43], [37,38,44,45], [39,40,46],
  [41,42], [43,44], [45,46]];

  stabXPRP5 := [
    SqStab(1,8), SqStab(3,10), SqStab(5,12), 
    SqStab(7,14), SqStab(9,16), SqStab(11,18), 
    SqStab(15,22), SqStab(19,25), 
    SqStab(21,27), SqStab(24,31), 
    SqStab(28,35), SqStab(30,37), SqStab(32,39), 
    SqStab(34,41), SqStab(36,43), SqStab(38,45)
  ];


stabPRP5 := [];
AddElementsOf(stabPRP5, MakeStab(stabZPRP5,"Z"));
AddElementsOf(stabPRP5, MakeStab(stabXPRP5,"X"));

# stab12x19RSC_1b_3VT3_1HT4 12x19 rotated surface code with single boundary and 4 twists, three vertical of length 3 a one horizontal of length 4

stabZ12x19RSC_1b_3VT3_1HT4 := [
  [1,2], [3,4], [5,6], [7,8], [9,10], [11,12],
  [1,13,14], SqStab(2,15), SqStab(4,17), SqStab(6,19), SqStab(8,21), SqStab(10,23), [12,25],
  [13,26], SqStab(16,28), SqStab(24,33),
  SqStab(26,35), SqStab(31,40), [34,43],
  [35,44], SqStab(37,47), SqStab(42,55),
  SqStab(44,57), SqStab(46,59), SqStab(48,61), SqStab(50,63), SqStab(52,65), SqStab(54,67), [56,69],
  [57,70], SqStab(58,71), SqStab(60,73),  SqStab(62,75), SqStab(64,77), SqStab(66,79), SqStab(68,81),
  SqStab(70,83), SqStab(72,85), SqStab(78,88), SqStab(80,90), [82,92],
  [83,93], SqStab(84,94), SqStab(86,96), SqStab(89,102), SqStab(91,104),
  [93,94,106], SqStab(95,107), SqStab(97,109), SqStab(99,111), SqStab(101,113), SqStab(103,115), [105,117],
  [106,107], [108,109], [110,111], [112,113], [114,115], [116,117]
];

stabX12x19RSC_1b_3VT3_1HT4 := [
  SqStab(1,14), SqStab(3,16), SqStab(5,18), SqStab(7,20), SqStab(9,22), SqStab(11,24),
  SqStab(13,26), SqStab(21,31),
  SqStab(28,37), SqStab(33,42),
  SqStab(35,44), SqStab(40,52), 
  SqStab(45,58), SqStab(47,60), SqStab(49,62), SqStab(51,64), SqStab(53,66), SqStab(55,68),
  SqStab(57,70), SqStab(59,72), SqStab(61,74), SqStab(63,76), SqStab(65,78), SqStab(67,80),
  SqStab(71,84), SqStab(73,86), SqStab(79,89), SqStab(81,91),
  SqStab(83,93), SqStab(85,95), SqStab(88,101), SqStab(90,103),
  SqStab(94,106), SqStab(96,108), SqStab(98,110), SqStab(100,112), SqStab(102, 114), SqStab(104,116) 
];

stabXZ12x19RSC_1b_3VT3_1HT4 := [
  [[15,16,28],[14,15,27]],
  [[27,36], [28,37]],
  [[37,46,47],[36,45,46]],
  [[17,18,29],[18,19,30]],
  [[30,39],[29,38]],
  [[38,48,49],[39,49,50]],
  [[19,20,30],[20,21,31]],
  [[31,40],[30,39]],
  [[39,50,51],[40,51,52]],
  [[23,24,33],[22,23,32]],
  [[32,41],[33,42]],
  [[42,54,55],[41,53,54]],
  [[87,97,98],[74,75,87]],
  [[75,75],[98,99]],
  [[99,100],[76,77]],
  [[77,78,88],[88,100,101]]
];

stab12x19RSC_1b_3VT3_1HT4 := [];
AddElementsOf(stab12x19RSC_1b_3VT3_1HT4, MakeStab(stabZ12x19RSC_1b_3VT3_1HT4,"Z"));
AddElementsOf(stab12x19RSC_1b_3VT3_1HT4, MakeStab(stabX12x19RSC_1b_3VT3_1HT4,"X"));
AddElementsOf(stab12x19RSC_1b_3VT3_1HT4, stabXZ12x19RSC_1b_3VT3_1HT4);


# Andrew's Example of a single twist with one boundary
Xstab:=[[1,4,5],[2,3,6,7],[7,8,11,12],[9,10,13,14],[15,16,20,21],[17,18,22],[19,20,23,24]];
Zstab:=[[1,4],[1,2,5,6],[2,3],[3,7,8],[4,5,9,10],[8,12],[9,13],[11,12,15,16],[13,14,17,18],[16,21],[17,22],[18,19,22,23],[23,24],[20,21,24]];


# 9 quibit surface code  - rotated with single twist
Xstab := [];

ST:=[
[1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0],
[0,1,1,0,1,0,0,0,1,1,0,1,0,0,0,0],
[0,0,0,1,0,1,1,0,0,0,0,0,1,0,0,1],
[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]] * One(GF(2));;

ST_heads:= [ 1, 2, 3, 4, 5, 6, 7, 0, 8, 9, 10, 11, 12, 0, 13, 14, 15, 0, 0, 0, 16, 0 ];

BS:=
[ [  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0 ],
  [  0,  1,  1,  0,  1,  0,  0,  0,  1,  1,  0,  1,  0,  0,  0,  0 ],
  [  0,  0,  0,  1,  0,  1,  1,  0,  0,  0,  0,  0,  1,  0,  0,  1 ],
  [  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 ],
  [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 ] ] *One(GF(2));

  # Drew's Diamond like example - but when center quibit is removed without decreasing the number of stabvilizers the groups becomes non-abelian

Xstab:=[[1,2], [3,4], [1,5], [2,3,6,7], [4,8,9], [5,6,10,11], [7,8,13], [9,14], [10,15], [11,16,17], [13,14,18,19], [15,16,20], [17,18,21,22], [19,23], [20,21], [22,23]];
Zstab:=[[1,2,5,6],[3,4,7,8],[6,7,11],[8,9,13,14],[10,11,15,16],[13,17,18],[16,17,20,21],[18,19,22,23]];

# Surface code - two boundaries 4 x 5 square with a single twist (not to noundary - two vertical quibits removed)
# Encoding [S1, S2, ... ,SN] = [[U1,V1], [U2,V2] ..., [UN,VN]] where UN and VN are indices (1... pos) of X and Z's respectively : e.g.
# [[1,2],[2,4]] = X1 Y2 Z4 - phase not included
stab45 := [[[],[1,6]],
[[],[2,3,7,8]],
[[],[4,5,9,10]],
[[],[6,7,11,12]],
[[],[10,14]],
[[],[11,15]],
[[],[13,14,17,18]],
[[],[15,16,19,20]],
[[],[18,23]],
[[],[19,24]],
[[],[20,21,25,26]],
[[],[22,23,27,28]],
[[1,2,6,7],[]],
[[3,2],[]],
[[3,4,8,9],[]],
[[4,5],[]],
[[9,10,13,14],[]],
[[11,12,15,16],[]],
[[17,18,22,23],[]],
[[19,20,24,25],[]],
[[25,26],[]],
[[21,22,26,27],[]],
[[27,28],[]],
[[7,8,12],[8,9,13]],
[[13,17],[12,16]],
[[16,20,21],[17,21,22]]];


# Surface code 
# Encoding [S1, S2, ... ,SN] = [[U1,V1], [U2,V2] ..., [UN,VN]] where UN and VN are indices (1... pos) of X and Z's respectively : e.g.
# [[1,2],[2,4]] = X1 Y2 Z4 - phase not included
stabRSC9 := [[[1,2],[]],
[[2,3,5,6],[]],
[[4,5,7,8],[]],
[[8,9],[]],
[[],[1,2,4,5]],
[[],[3,6]],
[[],[4,7]],
[[],[5,6,8,9]]];

stabRSSC11 := [[[1,2],[]],
[[10,11],[]],
[[4,6,7],[]],
[[2,3,4],[]],
[[8,9,10],[]],
[[5,6,8],[]],
[[],[1,5]],
[[],[7,11]],
[[],[3,4,7]],
[[],[2,4,6]],
[[],[6,8,10]],
[[],[5,8,9]]];

# Basic Tests (some)

# One: 
TestOne := function()
  local G_int, G;

  G_int := CreateGenMatrix(stabRSSC11,11);
  G := G_int*One(GF(2));
  return IsotropicHyperboicForm(G);

end;

# Two
TestTwo := function()
  local G_int, G;

  G_int := CreateGenMatrix(stabRSSC11,11);
  G := G_int*One(GF(2));

  return GaugeGroupCenter(G);
end;

# Three
# 
TestThree := function()
  local G_int, G;

  G_int := [[1,0,0,0],[0,1,1,0],[0,0,1,0],[1,0,1,0]];
  G := G_int*One(GF(2));

  return GaugeGroupCenter(G,true);
end;

# Four
# 
TestFour := function()
  local G_int, G;

  G_int := [[1,0,0,0],[0,1,1,0],[0,0,1,0],[1,0,1,0]];
  G := G_int*One(GF(2));

  return GaugeGroupCenter(G,false);
end;

# Five
TestFive := function()
  local G_int, G;

  G_int := CreateGenMatrix(stabRSSC11,11);
  G := G_int*One(GF(2));

  return IsotropicHyperboicBasis(IsotropicHyperboicForm(G));
end;


# Older Tests


TestStabPRP5 := function()  
  local gen_mat, G, qec_info;
  gen_mat := CreateGenMatrix(stabPRP5,46);
  G := gen_mat*One(GF(2));
  qec_info := rec(isstabilizergroup:=true);

  return GaugeGroupNormalizer(G, qec_info);
end;

Teststab12x19RSC_1b_3VT3_1HT4 := function()  
  local gen_mat, G, qec_info;
  gen_mat := CreateGenMatrix(stab12x19RSC_1b_3VT3_1HT4,117);
  G := gen_mat*One(GF(2));
  qec_info := rec(isstabilizergroup:=true);

  return GaugeGroupNormalizer(G, qec_info);
end;


#  Three
# Assumes that S is a Stabilizer group i.e. abelian
MakeRSC9Three := function()
  local S_int, S, S_ext, S_ext_hyb;

  S_int := CreateGenMatrix(stabRSC9,9);
  S := S_int*One(GF(2));

  S_ext := ExtendBasis(S);
  S_ext_hyb := GramSchmidt(S_ext);

  return RemoveHyperbolicBasisElements(S_ext_hyb, "HyperZ", [1..Size(S)]);
end;

# Four
MakeRSSC11Four := function()
  local G_int, G, G_ext, G_ext_hyp, G_tmp, distance_of_center, N_G, gauge_degrees, number_of_logical_qubits,Y,range, G_ext_com, orange, xrange, zrange;

  G_int := CreateGenMatrix(stabRSSC11,11);
  G := G_int*One(GF(2));

  G_tmp := GramSchmidt(G);
  # Note: technically distance of Z(G) is distance_of_center + 1 (as iI is a needed generator)
  distance_of_center := Size(G_tmp.center);
  gauge_degrees := Size(G_tmp.hyperX);

  G_ext := ExtendBasis(CombineComponents(G_tmp.hyperX, G_tmp.hyperZ, G_tmp.center));

  orange := [(2*Size(G_tmp.hyperX)+1)..Size(G_ext)];
  xrange := [1..Size(G_tmp.hyperX)];
  zrange := [Size(G_tmp.hyperX)+1..2*Size(G_tmp.hyperX)];

  G_ext_com := MakeSetsCommuteGeneral(G_ext, orange , G_ext, xrange, G_ext, zrange);
  G_ext_hyp := GramSchmidt(G_ext_com, 2*Size(G_tmp.hyperX)+1);

  range := [1+Size(G_tmp.hyperX)..distance_of_center+Size(G_tmp.hyperX)];
  N_G := RemoveHyperbolicBasisElements(G_ext_hyp, "HyperZ", range);
  number_of_logical_qubits := Size(N_G.hyperX)-gauge_degrees;

  return rec(num_qubits:=number_of_logical_qubits, center_distance:=distance_of_center, gauge_degrees:=gauge_degrees, normalizer:=N_G);
end;


MakeRSSC11Five := function()
  local G_int, G, center_size, gauge_size, G_enb, G_ext, orange, xrange, zrange, G_ext_com, G_hyp, range, N_G;

  G_int := CreateGenMatrix(stabRSSC11,11);
  G := G_int*One(GF(2));
  G := GramSchmidt(G);

  center_size := Size(G.center);
  gauge_size := Size(G.hyperX);

  G_enb := ExtendHyperbolicBasis(G);
  G_ext := ExtendBasis(CombineComponents(G_enb.hyperX, G_enb.hyperZ));

  orange := [2*Size(G_enb.hyperX)+1..Size(G_ext)];
  xrange := [1..Size(G_enb.hyperX)];
  zrange := [Size(G_enb.hyperX)+1..2*Size(G_enb.hyperX)];

  G_ext_com := MakeSetsCommuteGeneral(G_ext, orange, G_ext, xrange, G_ext, zrange);
  G_hyp := GramSchmidt(G_ext_com, 2*Size(G_enb.hyperX)+1);

  range := [gauge_size+1..gauge_size+center_size];
  N_G := RemoveHyperbolicBasisElements(G_hyp, "HyperZ", range);

  return rec(N_G:=N_G, gauge_degree:=gauge_size, logicals := Size(N_G.hyperX)-gauge_size);
end;


MakeRSSC11Six := function()
  local G_int, G, center_size, gauge_size, G_enb, G_ext, orange, xrange, zrange, G_ext_com, G_hyp, range,P_n;

  G_int := CreateGenMatrix(stabRSSC11,11);
  G := G_int*One(GF(2));
  G := GramSchmidt(G);

  center_size := Size(G.center);
  gauge_size := Size(G.hyperX);

  G_enb := ExtendHyperbolicBasis(G);
  G_ext := ExtendBasis(CombineComponents(G_enb.hyperX, G_enb.hyperZ));

  orange := [2*Size(G_enb.hyperX)+1..Size(G_ext)];
  xrange := [1..Size(G_enb.hyperX)];
  zrange := [Size(G_enb.hyperX)+1..2*Size(G_enb.hyperX)];

  G_ext_com := MakeSetsCommuteGeneral(G_ext, orange, G_ext, xrange, G_ext, zrange);
  P_n := GramSchmidt(G_ext_com, 2*Size(G_enb.hyperX)+1);

  return P_n;
end;
