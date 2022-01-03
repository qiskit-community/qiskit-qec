import numpy as np
import matplotlib.pyplot as plt
from pymatching import Matching
from plerco.pec_python.qec_code.rssc import RSSC

# Simulate RSSC with channel model

np.random.seed(1)

dscan = [3, 5, 7, 11, 15, 29]
pscan = [0.10, 0.12, 0.14, 0.16, 0.18]
samples = 10000


results = {k: [] for k in dscan}
for d in dscan:
    code = RSSC(d)
    Hz = np.array(code.z_gauge_matrix)
    Lz = np.array(code._to_matrix([code.logical_z]))
    m = Matching(Hz)
    plist = []
    for p in pscan:
        fail = 0
        noise_weights = []
        correction_weights = []
        for i in range(samples):
            noise = (np.random.rand(code.n) < p).astype(np.uint8)
            noise_weights.append(np.count_nonzero(noise))
            syn = Hz @ noise % 2
            c = m.decode(syn)
            correction_weights.append(np.count_nonzero(c))
            res = (c + noise) % 2
            syn0 = Hz @ res % 2
            if len(np.nonzero(syn0)[0]) > 0:
                raise Exception("error correction failed!")
            if (Lz @ res)[0] % 2 == 1:
                fail += 1
                if noise_weights[-1] <= (d - 1) / 2:
                    print(
                        "noise_weight = %d, correction_weight = %d"
                        % (noise_weights[-1], correction_weights[-1])
                    )
                    print(noise)
                    print(syn)
                    print(c)
                    print(res)
                    print(syn0)
                    print(Lz @ res)
        pf = float(fail) / float(samples)
        print(
            "%d: %f  -->  %f (%f, %f)"
            % (d, p, pf, np.mean(noise_weights), np.mean(correction_weights))
        )
        plist.append(pf)
    results[d] = plist
print(results)

dscan = [3, 5, 7, 11, 15, 29]

plt.figure(figsize=(12, 8))
plt.loglog(pscan, results[3], "bo-")
plt.loglog(pscan, results[5], "rx-")
plt.loglog(pscan, results[7], "gp-")
plt.loglog(pscan, results[11], "bo-")
plt.loglog(pscan, results[15], "rx-")
plt.loglog(pscan, results[29], "gp-")
plt.grid(which="major", axis="both")
plt.xlabel("p(fail)")
plt.ylabel("p")
# plt.show()
plt.savefig("channel-rssc-gauge.png")
