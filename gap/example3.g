#
# Example: 7-qubit code
#
Read("aut.g");

genmat := StringsToGF4(["xixixix", "ixxiixx", "iiixxxx", "ziziziz", "izziizz", "iiizzzz"]);
nq := 7;

t1 := NanosecondsSinceEpoch();
A := GF4AdditiveAutGroup(genmat);
t2 := NanosecondsSinceEpoch();
Print("time = ", Float(t2 - t1)/1.e9, " seconds\n");
Print("Order(A) = ", Order(A), "\n");
Print(reformat(A, nq), "\n");

Aperm := GF4AdditivePermAutGroup(genmat);
Print("Order(Aperm) = ", Order(Aperm), "\n");
Print(reformat(Aperm, nq), "\n");

Alocal := GF4AdditiveLocalAutGroup(genmat);
Print("Order(Alocal) = ", Order(Alocal), "\n");
Print(reformat(Alocal, nq), "\n");
