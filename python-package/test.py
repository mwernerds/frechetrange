import numpy as np;
import frechet;


if __name__=="__main__":
    print("Create some trajectory");
    dec = frechet.FrechetDecider()
    t1 = np.array([[1,0],[2,0],[4,0],[5,0]])
    t2 = np.array([[1,0],[3,3],[5,0]])
    print(t1)
    print(t2)

    print("DV @ 2.99: %d" % (dec.decide_dv(t1,t2,2.99)))
    print("DV @ 3.01: %d" % (dec.decide_dv(t1,t2,3.01)))
    print("Decide @ 3.01: %d" % (dec.decide(t1,t2,3.01,"duetschvahrenhold")))
    print("Decide @ 3.01: %d" % (dec.decide(t1,t2,3.01,"baldusbringman")))

    

    
