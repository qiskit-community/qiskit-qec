from plerco.pec_python.simulation_runner.runner import Runner

pin = [0.0008, 0.001, 0.003]
r = Runner()
for code in ["hhc"]:
    for d in [3, 5, 7, 9, 11]:
        for rounds in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            r.add_experiment(code, "slow", d, rounds, "zzxx", "z", pin)
            r.add_experiment(code, "slow", d, rounds, "zzxx", "x", pin)
r.run_experiments()
print("DONE")
print(r.results)
