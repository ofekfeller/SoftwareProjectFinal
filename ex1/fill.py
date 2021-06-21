


def find_sub(S):
    P = []
    i=0
    sum= sum(S)/2
    return rec(S,P,i,sum)



def rec(S, P , index, sum):
    if (index==len(S)): return 0
    #blablabla
