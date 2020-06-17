import numpy as np;
import frechet

# this is python glue code to clean up towards our format
def frechet_distance(t1,t2):
    # safely move to a python list of list with float values inside from any sufficient iteratable
    lt1 = [list(map(float,list(x))) for x in t1]
    lt2 = [list(map(float,list(x))) for x in t2]
    print(lt1)
    print(lt2)
    return frechet.frechet_distance(lt1,lt2)

if __name__=="__main__":
    # Integer lists, not supported by the extension, but through above glue code
    t1 = [[x,y+1] for x,y in zip(range(10),range(10))]
    t2 = [[x+1,y+2] for x,y in zip(range(10),range(10))]
    print("Integer list FD: %f" % (frechet_distance(t1,t2)))
    
    # Numpy arrays as well
    nt1 =np.array(t1)
    nt2 =np.array(t2)
    print("Numpy int FD: %f" % (frechet_distance(nt1,nt2)))

    nt1 =np.array(t1, dtype=np.double)
    nt2 =np.array(t2,dtype=np.double)
    print("Numpy float FD: %f" % (frechet_distance(nt1,nt2)))

    # And at least a single test (proof by example ;-)
    t1 = np.array([[1,0],[2,0],[4,0],[5,0]], np.double)
    t2 = np.array([[1,0],[3,3],[5,0]])
    # Distance is three (draw the triangle)
    print("Our implementation knows that 3=%f" % (frechet_distance(t1,t2)))

    
