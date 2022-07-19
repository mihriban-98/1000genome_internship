
def adict(dict1,dict2):
    for i in dict1.keys():
        for j in dict2.keys():
            if i == j:
                dict1[i] = dict1[i] + dict2[j]

    for j in dict2.keys():
        #print(dict1)
        if j not in dict1.keys():
            a = {j : dict2[j]}
            dict1.update(a)
    return dict1


