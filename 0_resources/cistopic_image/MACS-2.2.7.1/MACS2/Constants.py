MACS_VERSION = "2.2.7.1"
FILTERDUP_VERSION = "1.0.0 20140616"
RANDSAMPLE_VERSION = "1.0.0 20120703"
MAX_PAIRNUM = 1000
MAX_LAMBDA  = 100000
FESTEP      = 20
BUFFER_SIZE = 100000                   # np array will increase at step of 1 million items

from array import array

if array('h',[1]).itemsize == 2:
    BYTE2 = 'h'
else:
    raise Exception("BYTE2 type cannot be determined!")

if array('H',[1]).itemsize == 2:
    UBYTE2 = 'H'
else:
    raise Exception("UBYTE2 (unsigned short) type cannot be determined!")

if array('i',[1]).itemsize == 4:
    BYTE4 = 'i'
elif array('l',[1]).itemsize == 4:
    BYTE4 = 'l'
else:
    raise Exception("BYTE4 type cannot be determined!")

if array('f',[1]).itemsize == 4:
    FBYTE4 = 'f'
elif array('d',[1]).itemsize == 4:
    FBYTE4 = 'd'
else:
    raise Exception("FBYTE4 type cannot be determined!")
