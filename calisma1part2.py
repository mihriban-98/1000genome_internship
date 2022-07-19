import pandas as pd
import numpy as np
from deneme import adict

#Needed names
continents = ["AFR","AMR","EAS","EUR","SAS"]
rslerlist = [ "rs1815739", "rs4644994", "rs4253778", "rs1042713", "rs8192678", "rs2070744", "rs11549465", "rs1800012", "rs12722", "rs1049434",
        "rs1800795", "rs6265", "rs4680"]

#Dataframes created for different ethnicities
AFR = pd.DataFrame()
AMR = pd.DataFrame()
EAS = pd.DataFrame()
EUR = pd.DataFrame()
SAS = pd.DataFrame()

#Dataframes filled with according data
for i in rslerlist:
    for j in continents:
        a = pd.read_csv("data/ethnicity/{}/{}.csv".format(i,j),names=rslerlist)
        a = a.iloc[:,1]
        #print(a)
        if j == "AFR":
            AFR = pd.concat([AFR,a],axis=1)
            #print(AFR)
        elif j == "AMR":
            AMR = pd.concat([AMR,a], axis=1)
        elif j == "EAS":
            EAS = pd.concat([EAS, a], axis=1)
        elif j == "EUR":
            EUR = pd.concat([EUR,a], axis=1)
        elif j == "SAS":
            SAS = pd.concat([SAS,a], axis=1)

#Column names changed to correct names
AFR.columns = rslerlist
AMR.columns = rslerlist
EAS.columns = rslerlist
EUR.columns = rslerlist
SAS.columns = rslerlist

#print(AFR.rs1815739)

#print(AFR.isna().sum().sum())
#Frquency calculator for homo, hetero and genotype genes

total_number = []
for df in [AFR, AMR, EAS, EUR, SAS]:
    total_number.append(df.count(axis=1))

AFR_tot = (int(len(total_number[0])) * 13)
AMR_tot = (int(len(total_number[1])) * 13)
EAS_tot = (int(len(total_number[2])) * 13)
EUR_tot = (int(len(total_number[3])) * 13)
SAS_tot = (int(len(total_number[4])) * 13)


def allel_freq_calc():

    AFR_frequency = {}
    AMR_frequency = {}
    EAS_frequency = {}
    EUR_frequency = {}
    SAS_frequency = {}

    AFR_values = AFR.rs1815739.value_counts().to_dict()
    AMR_values = AMR.rs1815739.value_counts().to_dict()
    EAS_values = EAS.rs1815739.value_counts().to_dict()
    EUR_values = EUR.rs1815739.value_counts().to_dict()
    SAS_values = SAS.rs1815739.value_counts().to_dict()

    values_list = [AFR_values, AMR_values, EAS_values, EUR_values, SAS_values]

    #Number of different genotypes on different ethnicities
    for i in AFR.columns[1:14]:

        AFR_values = adict(AFR_values, AFR[i].value_counts().to_dict())
        AMR_values = adict(AMR_values, AMR[i].value_counts().to_dict())
        EAS_values = adict(EAS_values, EAS[i].value_counts().to_dict())
        EUR_values = adict(EUR_values, EUR[i].value_counts().to_dict())
        SAS_values = adict(SAS_values, SAS[i].value_counts().to_dict())

    #print(AFR_values)

    heterovals = [["A|C","C|A"],["A|T","T|A"],["A|G","G|A"],["C|T","T|C"],["G|C","C|G"],["T|G","G|T"]]

    for i in values_list:

        for j,k in heterovals:
            if j in i.keys():
                if k in i.keys():
                    #print(i[k])
                    new_val = i[j] + i[k]
                    del i[j],i[k]
                    new_row = {"{}".format(j): new_val }
                    i.update(new_row)

        #print(i)

    print(AFR_values)

    for i in values_list:

        for j in ["A", "C", "T", "G"]:
            for k in ["A", "C", "T", "G"]:
                if "{}|{}".format(j,k) in i.keys():

                    if i == AFR_values:

                        new_vals = {"{}|{}".format(j, k): ((i["{}|{}".format(j, k)] / AFR_tot) * 100)}
                        AFR_frequency.update(new_vals)
                    elif i == AMR_values:
                        new_vals = {"{}|{}".format(j, k): ((i["{}|{}".format(j, k)] / AMR_tot) * 100)}
                        AMR_frequency.update(new_vals)
                    elif i == EAS_frequency:
                        new_vals = {"{}|{}".format(j, k): ((i["{}|{}".format(j, k)] / EAS_tot) * 100)}
                        EAS_frequency.update(new_vals)
                    elif i == EUR_values:
                        new_vals = {"{}|{}".format(j, k): ((i["{}|{}".format(j, k)] / EUR_tot) * 100)}
                        EUR_frequency.update(new_vals)
                    elif i == SAS_values:
                        new_vals = {"{}|{}".format(j, k): ((i["{}|{}".format(j, k)] / SAS_tot) * 100)}
                        SAS_frequency.update(new_vals)

    return AFR_frequency, AMR_frequency, EAS_frequency, EUR_frequency, SAS_frequency


AFR_frequency, AMR_frequency, EAS_frequency, EUR_frequency, SAS_frequency = allel_freq_calc()

print(AFR_frequency)








