import re
import time
import timeit
import gmpy2
import math
import multiprocessing as mp
import GaussElimination
import sqrt_modulo

__author__ = 'son'

gmpy2.get_context().precision=500*8


n50b = "543971971236630195076898180273"
n60b = "459344432904829185748347139878533089"
n70b = "810814527278919159798277223576982175301717"
n80b = "662875573009490731045345994740327457245824480677"
n25b = 767971435180333
test = 90283
N = gmpy2.mpz(n70b)
n_smooth = 1500
offset = 800000000
start = 0
listN = []

start = gmpy2.mpz(gmpy2.ceil(gmpy2.sqrt(N)))

def init():
    for j in range(0, offset):
        listN.append(gmpy2.sub(gmpy2.mul(gmpy2.add(start, j), gmpy2.add(start, j)), N))
    return start


def read_prime():
    list_primes = []
    with open("primes1000.txt") as fp:
        for line in fp:
            list_primes.extend([int(s) for s in re.split(r'\s+', line.strip())])
    print("list prime:")
    print(list_primes)
    return list_primes


def read_10k_prime():
    list_primes = []
    with open("primes10000.txt") as fp:
        for line in fp:
            list_primes.extend([int(s) for s in re.split(r'\s+', line.strip())[1:]])

    legendre_prime = [2]
    for i in range(1, len(list_primes)):
        if gmpy2.legendre(N, list_primes[i]) == 1:
            legendre_prime.append(list_primes[i])

    print("len list prime:" + str(len(legendre_prime)))
    print(legendre_prime[:20])
    return legendre_prime


def read_100k_prime():
    list_primes = []
    with open("100k_primes.txt") as fp:
        for line in fp:
            list_primes.append(int(line.strip()))

    legendre_prime = [2]
    for i in range(1, len(list_primes)):
        if gmpy2.legendre(N, list_primes[i]) == 1:
            legendre_prime.append(list_primes[i])

    print("len list prime:" + str(len(legendre_prime)))
    print(legendre_prime[:20])
    return legendre_prime


def count_nth_smooth(listNN):
    count = 0
    for n in listNN:
        if n == 1:
            count += 1
    return count


def create_primes_index(list_primes):
    list_start_index = []
    for i in range(0, n_smooth + 1):
        p = list_primes[i]
        sqr1 = sqrt_modulo.sqrt_modular(N, p)
        sqr2 = p - sqr1
        start_index1 = gmpy2.f_mod(gmpy2.sub(sqr1 - start), p)
        start_index2 = gmpy2.f_mod(gmpy2.sub(sqr2 - start), p)
        list_start_index.append([start_index1, start_index2])
    return list_start_index


def get_smoothness(_start, delta, list_primes, pos, output):
    print("start process " + str(pos))

    segment = math.ceil(delta / 10)
    V = []
    for k in range(0, 10):
        start_p = pos * delta + k * segment
        end_p = pos * delta + (k + 1) * segment
        if end_p > offset:
            end_p = offset

        listNN = []
        for j in range(start_p, end_p):
            listNN.append(gmpy2.sub(gmpy2.mul(gmpy2.add(_start, j), gmpy2.add(_start, j)), N))

        for i in range(0, n_smooth + 1):
            p = list_primes[i]
            index = 0
            for index in range(0, len(listNN)):
                if gmpy2.f_mod(listNN[index], p) == 0:
                    break
            index += start_p

            start_index1 = gmpy2.f_mod(index, p)

            for n in range(index - start_p, len(listNN), p):
                while gmpy2.f_mod(listNN[n], p) == 0:
                    listNN[n] = gmpy2.div(listNN[n], p)

            start_index2 = gmpy2.f_mod(gmpy2.sub(gmpy2.sub(p, gmpy2.f_mod(gmpy2.add(_start, start_index1), p)), _start),p)
            index2 = gmpy2.sub(p, gmpy2.f_mod(gmpy2.sub(start_p, start_index2), p))

            if start_index2 != start_index1:
                for n in range(index2, len(listNN), p):
                    while gmpy2.f_mod(listNN[n], p) == 0:
                        listNN[n] = gmpy2.div(listNN[n], p)

        for i in range(0, len(listNN)):
            if listNN[i] == 1:
                start1 = gmpy2.add(_start, start_p)
                num = gmpy2.sub(gmpy2.mul(gmpy2.add(start1, i), gmpy2.add(start1, i)), N)
                V.append(num)

    output.put((pos, V))
    print("end process " + str(pos))
    # print(listNN)


def sieve_of_eratosthenes(_start, list_primes):
    print("tap cac so phai xet:"+str(offset))
    # print(listN[:20])
    print("n_smooth: "+str(n_smooth))
    delta = math.ceil(offset / 4)

    output = mp.Queue()
    processes = [mp.Process(target=get_smoothness, args=(_start, delta, list_primes, x, output)) for x in range(0, 4)]
    for p in processes:
        p.start()
    results = [output.get() for p in processes]
    results.sort()
    tem = []
    for r in results:
        tem = tem + r[1]
    results = tem
    # print("result")
    # print(results)
    print("num row:")
    print(len(results))
    return results


def create_gauss_matrix(V, list_primes):
    A = [[]]
    for i in range(0, len(V)):
        num = V[i]
        vector = []
        for p in range(0, n_smooth + 1):
            value = 0
            while gmpy2.c_mod(num, list_primes[p]) == 0:
                num = gmpy2.div(num, list_primes[p])
                value ^= 1
            vector.append(value)
        A.append(vector)

    print("cac so co the tao tich binh phuong: " + str(len(V)))
    print(V[:20])
    return A[1:]


# phan tich gauss
def gauss_elimination(A):
    print("ma tran cho phan tich gauss:" + str(len(A)) + "  " + str(len(A[0])))
    # tao ma tran chuyen vi
    B = []
    for j in range(0, len(A[0])):
        v = [0]*len(A)
        for i in range(0,len(A)):
            v[i] = A[i][j]
        B.append(v)
    print(A[:10])
    print("ma tran chuyen vi: ")
    print(B[:10])
    A = []
    # phan tich gauss
    r = GaussElimination.new_sort(B)
    return r


def get_a_result(C, V, n_smooths, list_primes):
    row = len(C)
    col = len(C[0])
    result = [-1]*col

    for i in range(0, col):
        if C[i][i] == 0:
            C = C[:i] + [[0]*col] + C[i:]

    first = 1
    for i in range(col - 1, -1, -1):
        if C[i][i] == 0:
            if first == 1:
                result[i] = 1
                first = 0
            else:
                result[i] = 0
        else:
            result[i] = 0
            for j in range(i + 1, col):
                result[i] ^= C[i][j]*result[j]
    print("result:")
    print(result)

    sqr = gmpy2.mpz(1)
    product = gmpy2.mpz(1)
    ex = [0]*(n_smooths + 1)

    for i in range(0,  len(V)):

        if result[i] == 0:
            continue
        num = V[i]
        for j in range(0, n_smooths + 1):
            n_prime = list_primes[j]
            while gmpy2.f_mod(num, n_prime) == 0:
                ex[j] += 1
                num = gmpy2.div(num, n_prime)
        product = gmpy2.f_mod(gmpy2.mpz(gmpy2.mul(product, gmpy2.sqrt(gmpy2.add(V[i], N)))), N)
    print("ex")
    print(ex)
    # tinh square
    for i in range(0, n_smooths + 1):
        num = list_primes[i]
        e = ex[i] / 2
        sqr = gmpy2.f_mod(gmpy2.mul(sqr, gmpy2.powmod(gmpy2.mpz(num), gmpy2.mpz(e), N)), N)

    print("product square")
    print(sqr)
    print("product")
    print(product)

    # print((gmpy2.sub(product, sqr)))
    # print(gmpy2.mpz(gmpy2.sub(product, sqr)))

    p = gmpy2.gcd(gmpy2.mpz(gmpy2.sub(product, sqr)), N)
    q = gmpy2.gcd(gmpy2.mpz(gmpy2.add(product, sqr)), N)
    print("p, q")
    print(str(p) + " is prime: " + str(gmpy2.is_prime(p)))
    print(str(q) + " is prime: " + str(gmpy2.is_prime(q)))

    return result


def get_b_smooth():
    t = gmpy2.sqrt(gmpy2.mul(gmpy2.log(N), gmpy2.log(gmpy2.log(N))))
    b= gmpy2.exp(gmpy2.mul(0.5, t))
    print("bsmooth advice: " + str(b))


def main():
    start_time = timeit.default_timer()
    # init(start)
    print("start")
    print(start)

    get_b_smooth()
    list_primes = read_100k_prime()
    V = sieve_of_eratosthenes(start, list_primes)

    gauss_matrix = create_gauss_matrix(V, list_primes)
    result = gauss_elimination(gauss_matrix)
    get_a_result(result, V, n_smooth, list_primes)

    end_time = timeit.default_timer()
    print(end_time - start_time)

if __name__ == '__main__':
    main()
