import numpy


def mean(vector):
    return float(sum(vector))/float(len(vector))


# Because of problems on same architectures use of unsigned integers is avoided.
def derive_dtype(n):
    """ Derive datatype for storing 0-based rankings for a given set length. """
    if n <= 2**15:
        # Range int16: -2^15 (= -32768) to 2^15 - 1 (= 32767).
        return numpy.int16
    else:
        # Range int32: -2^31 (= -2147483648) to 2^31 - 1 (= 2147483647).
        return numpy.int32
