def pearson(s1,s2):
    """ Computes the Pearson's correlation coeffient between s1 and s2 """
    d = len(s1)
    ki = sum(s1)
    kj = sum(s2)
    expR = (ki*kj)/d 
    pi = ki/d
    pj = kj/d
    upsiloni = sqrt(pi*(1-pi))
    upsilonj = sqrt(pj*(1-pj))
    r = sum((s1+s2)==2)
    similarity = (r-expR)/(d*upsiloni*upsilonj)
    if (ki==d and kj==d) or (ki==0 and kj==0):
        similarity = 1
    elif isnan(similarity):
        similarity = 0
   return(similarity)

def stab_index(A,pearson):
    """ Computes the average pairwise similarities between the rows of A """
    M = size(A,1)
    stability = 0
    for i in range(len(M)):
        for j in range(len(M)):
            if i != j:
                stability = stability + sim(A[i], A[j:])
    stability = stability/(M*(M-1))
    return(stability)
