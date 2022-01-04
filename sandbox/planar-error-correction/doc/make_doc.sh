#!/bin/bash
modules="css_code paulinoisemodel epselector errorpropagator faultenumerator hhc hhc_circuit circuit_matching_decoder hhc_decoder indexer rssc rssc_circuit rssc_decoder config"
for mod in ${modules}
do
pdoc --output-dir=doc --force ${mod}
done
