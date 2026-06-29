from bitarray import bitarray
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def hamming_distance(bits1, bits2):
    cdef Py_ssize_t n = len(bits1)
    assert n == len(bits2)
    cdef Py_ssize_t i
    cdef int distance = 0
    for i in range(0, n, 2):
        if bits1[i] != bits2[i] or bits1[i+1] != bits2[i+1]:
            distance += 1
    return distance 