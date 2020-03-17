

# this code intends to test the Ohm' module




import hecklerun as hr
import ohm



def ohmtest():
    r = hr.HeckleRun('r','/Volumes/drobo/nico/asymmetric/062/')
    ol = ohm.Ohm(r, 50)

    return [r,ol]




def main():
    "super main()"
    ohmtest()

if __name__ == '__main__':
    main()





