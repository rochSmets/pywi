#!/opt/local/bin/python2.7
# encoding: utf-8

import sys


def GetMemory(nx,nz,np):
    return 4*(10*(nx*nz) + 1*(nx+1)*(nz+1) + 5*np +11*nx*nz*2 + 2*nx*nz*2*4)




def main():
    mem =  GetMemory(int(sys.argv[1]),
                    int(sys.argv[2]),
                    int(sys.argv[3]))

    memKB = mem/1024.
    memMB = memKB/1024.
    memGB = memMB/1024.


    print 'memory usage (KB,MB,GB) : ',(memKB,memMB,memGB)



if __name__ == '__main__':
    main()


