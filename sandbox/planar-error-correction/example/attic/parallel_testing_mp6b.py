# process class in multiprocessing

import time
import multiprocessing


def multiprocessing_func(n):
    num = 1

    while n >= 1:
        num *= n
        n = n - 1

    print(num)


if __name__ == "__main__":
    starttime = time.time()
    processes = []
    for i in range(0, 1000):
        p = multiprocessing.Process(target=multiprocessing_func, args=(i,))
        processes.append(p)
        p.start()

    for process in processes:
        process.join()

    print("That took {} seconds".format(time.time() - starttime))
