# pool class in multiprocessing

import time
import multiprocessing as mp


def multiprocessing_func(n):
    num = 1

    while n >= 1:
        num *= n
        n = n - 1

    print(num)


if __name__ == "__main__":

    starttime = time.time()
    pool = mp.Pool()
    pool.map(multiprocessing_func, range(0, 1000))
    pool.close()
    print("That took {} seconds".format(time.time() - starttime))
