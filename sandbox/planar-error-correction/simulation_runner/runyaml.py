"""Script to run simulations from a YAML file."""
import time
import logging
import argparse
import signal
import os
import math
from datetime import datetime


from qiskit import execute, Aer
from qiskit.providers.aer import noise
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import *


from plerco.pec_python.qec_code.rssc import RSSC

from plerco.pec_python.qec_circuit.rssc_circuit import RSSCCircuit
from plerco.pec_python.qec_decoder.rssc_decoder import RSSCDecoder

from plerco.pec_python.qec_code.hhc import HHC
from plerco.pec_python.qec_circuit.hhc_circuit import HHCCircuit
from plerco.pec_python.qec_decoder.hhc_decoder import HHCDecoder

from plerco.pec_python.qec_noise_model.paulinoisemodel import PauliNoiseModel

from plerco.pec_python.configuration.config import config_load_all


exit_now = False  # flag to let us exit cleanly


def signal_handler(sig, frame):
    """SIGINT handler to let us exit cleanly."""
    global exit_now
    exit_now = True


def required_key(d, key):
    """Test if key is in dict d.

    Raise exception if not found.
    """
    if key not in d:
        raise Exception("key '%s' not found" % key)


def get_config(config):
    """Test and return values from a config."""
    # info section
    info = config["info"]
    # code section
    if config["code"] is None:
        raise Exception("missing 'code' section")
    required_key(config["code"], "distance")
    d = config["code"]["distance"]
    required_key(config["code"], "family")
    family = config["code"]["family"]
    # circuit section
    if config["circuit"] is None:
        raise Exception("missing 'circuit' section")
    required_key(config["circuit"], "rounds")
    rounds = config["circuit"]["rounds"]
    required_key(config["circuit"], "basis")
    basis = config["circuit"]["basis"]
    required_key(config["circuit"], "initial_state")
    initial_state = config["circuit"]["initial_state"]
    # simulation section
    if config["simulator"] is None:
        raise Exception("missing 'simulator' section")
    required_key(config["simulator"], "shots_increment")
    shots_increment = config["simulator"]["shots_increment"]
    required_key(config["simulator"], "shots_max")
    shots_max = config["simulator"]["shots_max"]
    required_key(config["simulator"], "stopping_threshold")
    stopping_threshold = config["simulator"]["stopping_threshold"]
    required_key(config["simulator"], "seed")
    seed = config["simulator"]["seed"]
    pin = None  # we don't require pin
    if "pin" in config["simulator"]:
        pin = config["simulator"]["pin"]
    return (
        info,
        d,
        family,
        rounds,
        basis,
        initial_state,
        shots_increment,
        shots_max,
        stopping_threshold,
        seed,
        pin,
    )


def sim_one(
    circ,
    model,
    decoder,
    basis,
    initial_state,
    seed,
    shots_increment,
    shots_max,
    stopping_threshold,
):
    """Simulate one circuit at error probably p.

    circ = quantum circuit
    model = pauli noise model
    decoder = decoder object
    basis = 'x' or 'z' prep/meas basis
    initial_state = '+' or '-' eigenvalue
    seed = simulator seed value
    shots_increment = number of shots per simulator call
    shots_max = largest number of shots to take
    stopping_threshold = number of failures before stopping

    Return (shots, pfail, sigma)
    """
    decoder.update_edge_weights(model)
    # Run chunks of shots_increment until we see enough failures
    total_shots = 0
    total_failures = 0
    while (
        total_failures < stopping_threshold and total_shots < shots_max and not exit_now
    ):
        tic = time.perf_counter()
        counts_n = run_circuit(circ, shots_increment, model, seed)
        toc = time.perf_counter()
        print(
            f"simulate %d circuits: {toc - tic:0.4f} seconds" % shots_increment,
            flush=True,
        )
        tic = time.perf_counter()
        failures = sum(
            map(
                lambda x: decode_outcome(x[0], x[1], decoder, basis, initial_state),
                counts_n.items(),
            )
        )
        toc = time.perf_counter()
        print(f"decode outcomes: {toc - tic:0.4f} seconds", flush=True)
        total_failures += failures
        total_shots += shots_increment
        seed += 1000
        print("   --> %d / %d failed" % (total_failures, total_shots), flush=True)
    pfail = float(total_failures) / float(total_shots)
    sig = sigma(pfail, total_shots)
    print("* pfail = %f +- %f (%d shots)" % (pfail, sig, total_shots), flush=True)
    return total_shots, pfail, sig


def run_circuit(qc, shots, model, seed):
    """Simulate a circuit with the stabilizer simulator.

    qc = QuantumCircuit to simulate
    shots = number of shots to run
    model = noise model to use
    seed = seed value for random number generator

    Return the dictionary of counts.
    """
    qasmsim = Aer.get_backend("qasm_simulator")
    if model is None:
        result = execute(
            qc,
            qasmsim,
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
    else:
        aer_model = model.as_aer_noise_model()
        result = execute(
            qc,
            qasmsim,
            noise_model=aer_model,
            basis_gates=aer_model.basis_gates,
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
    return result.get_counts(qc)


def sigma(p, shots):
    """Compute the sample standard deviation.

    p = Bernoulli parameter
    shots = number of samples

    Return the value.
    """
    return math.sqrt(p * (1 - p) / float(shots - 1))


def decode_outcome(outcome, num, decoder, basis, initial_state):
    """Decode an outcome to decide failure.

    outcome = outcome string from simulator
    num = number of occurrences of this outcome
    decoder = decoder object
    basis = prep/meas basis
    initial_state = initial eigenstate

    Return number of failures associated to this outcome.
    """
    # convert to list of integers
    reversed_outcome = list(map(int, outcome[::-1]))
    corrected_outcomes = decoder.process(reversed_outcome)
    if basis == "z":
        fail = decoder.code.logical_x_error(corrected_outcomes)
    elif basis == "x":
        fail = decoder.code.logical_z_error(corrected_outcomes)
    if initial_state == "-":
        fail[0] = (fail[0] + 1) % 2
    return num if bool(fail[0]) else 0


# Parse the command line arguments
parser = argparse.ArgumentParser(
    prog="runyaml",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Run a set of simulations.",
)
parser.add_argument("--log", type=str, default=None, help="debug logging to file LOG")
parser.add_argument(
    "--prefix", type=str, default="", help="filename prefix for output files"
)
parser.add_argument("yaml_file", type=str, help="input YAML 1.0 file")
args = parser.parse_args()

# Start logging
if args.log is not None:
    logging.basicConfig(filename=args.log, filemode="w", level=logging.DEBUG)

# Load YAML files
configlist = config_load_all(args.yaml_file)

prefix = args.prefix
if prefix != "":
    print("prefix = %s" % prefix)
cur_path = os.path.dirname(__file__)
yaml_base = os.path.basename(args.yaml_file)
yaml_head = os.path.splitext(yaml_base)[0]

# Set signal handler to give us a chance to clean up
signal.signal(signal.SIGINT, signal_handler)

for config in configlist:

    (
        info,
        d,
        family,
        rounds,
        basis,
        initial_state,
        shots_increment,
        shots_max,
        stopping_threshold,
        seed,
        pin,
    ) = get_config(config)
    assert family in ["RSSC", "HHC"]
    model = PauliNoiseModel(config["noise"])

    # Set up permitted iterables
    d_values = d if type(d) is list else [d]
    if pin is None:
        p_values = None  # nothing to iterate
    else:
        p_values = pin if type(pin) is list else [pin]
    rounds_values = rounds if type(rounds) is list else [rounds]

    for d_iter in d_values:
        for rounds_iter in rounds_values:

            print("description: %s" % info)
            print("distance %d, rounds %d" % (d_iter, rounds_iter))
            print("constructing circuits and decoders ... ", flush=True, end="")
            config["circuit"]["rounds"] = rounds_iter
            tic = time.perf_counter()
            if family == "RSSC":
                code = RSSC(d_iter)
                gen = RSSCCircuit(code, config)
                circ = gen.syndrome_measurement()
                decoder = RSSCDecoder(code, circ, model, config)
            elif family == "HHC":
                code = HHC(d_iter)
                gen = HHCCircuit(code, config)
                circ = gen.syndrome_measurement()
                decoder = HHCDecoder(code, circ, model, config)
            toc = time.perf_counter()
            print(f"{toc - tic:0.4f} seconds", flush=True)

            result = {}
            if p_values is None:
                # Nothing to iterate, run once
                print("-- running", flush=True)
                shots, pf, s = sim_one(
                    circ,
                    model,
                    decoder,
                    basis,
                    initial_state,
                    seed,
                    shots_increment,
                    shots_max,
                    stopping_threshold,
                )
                result["shots"] = shots
                result["pf"] = pf
                result["sigma"] = s
                curtime = datetime.now().strftime("%d-%m-%Y-%H:%M:%S")
                result["timestamp"] = curtime
            else:
                for p in p_values:
                    print("-- running %f" % p, flush=True)
                    model.set_scaled_error_probabilities(p)
                    shots, pf, s = sim_one(
                        circ,
                        model,
                        decoder,
                        basis,
                        initial_state,
                        seed,
                        shots_increment,
                        shots_max,
                        stopping_threshold,
                    )
                    curtime = datetime.now().strftime("%d-%m-%Y-%H:%M:%S")
                    result[p] = {
                        "shots": shots,
                        "pf": pf,
                        "sigma": s,
                        "timestamp": curtime,
                    }
                    if exit_now:
                        break

            # Store result into configuration
            config["result"] = result

            curtime = datetime.now().strftime("%d-%m-%Y-%H:%M:%S")
            filename = "%s%s-%s.yaml" % (prefix, yaml_head, curtime)
            filepath = os.path.relpath(filename, cur_path)
            config.save(filepath)

            if exit_now:
                break
        if exit_now:
            break
    if exit_now:
        break

if exit_now:
    print("----- caught SIGINT -----")

print("DONE")
