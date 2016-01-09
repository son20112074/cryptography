__author__ = 'son'


B = [[1, 1, 0, 0, 1, 1],
     [0, 1, 0, 0, 1, 1],
     [1, 0, 1, 0, 0, 0],
     [0, 0, 1, 0, 0, 0],
     [0, 0, 0, 0, 0, 0],
     [0, 0, 0, 1, 0, 1]]


# sắp xếp:
def sort(A):

    m = len(A)
    n = len(A[0])
    k = 0
    for i in range(0, n):
        # find hang 1 thu k dau tien
        if i % 100 == 0:
            print(i)
        j = 0
        for j in range(k, m):
            if A[j][i] == 1:
                break

        if A[j][i] == 0:
            continue

        # swap hang k va hang j
        tem = A[k]
        A[k] = A[j]
        A[j] = tem

        for j in range(k+1, m):
            if A[j][i] == 1:
                A[j] = [x ^ y for x, y in zip(A[j], A[k])]

        k += 1
    return A


def new_sort(a):
    m = len(a)
    n = len(a[0])

    list_matrix = []
    for i in range(0, len(a)):
        l = []
        for j in range(0, len(a[i])):
            if a[i][j] == 1:
                l.append(j)
        list_matrix.append(l)
    a = []
    k = 0
    print("phan tich")
    for i in range(0, n):
        if i % 100 == 0:
            print(i)
        # find hang 1 thu k dau tien
        j = 0
        for j in range(k, m):
            if not list_matrix[j]:
                continue
            if list_matrix[j][0] == i:
                break
        if not list_matrix[j]:
                continue
        if list_matrix[j][0] != i:
            continue

        # swap hang k va hang j
        tem = list_matrix[k]
        list_matrix[k] = list_matrix[j]
        list_matrix[j] = tem

        for j in range(k+1, m):
            # xor hang k va hang j
            if not list_matrix[j]:
                continue
            if list_matrix[j][0] == i:
                tem = list_matrix[j]
                for t in list_matrix[k]:
                    try:
                        list_matrix[j].index(t)
                        tem.remove(t)
                    except ValueError:
                        tem.append(t)
                tem.sort()
                list_matrix[j] = tem

        k += 1
    print("xong")
    a = [[0]*n]*m
    for t in range(0, len(list_matrix), 1):
        tem = [0]*n
        for i in list_matrix[t]:
            tem[i] = 1
        a[t] = tem
    list_matrix = []
    return a

if __name__ == "__main__":
    print(new_sort(B))
    print(sort(B))

