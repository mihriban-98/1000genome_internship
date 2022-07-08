import pandas as pd
import numpy as np
import os
import csv

rs1815739 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs1815739.csv")
rs4644994 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs4644994.csv")
rs4253778 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs4253778.csv")
rs1042713 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs1042713.csv")
rs8192678 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8192678.csv")
rs2070744 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs2070744.csv")
rs11549465 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs11549465.csv")
rs1800012 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs1800012.csv")
rs12722 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs12722.csv")
rs1049434 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs1049434.csv")
rs1800795 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs1800795.csv")
rs6265 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs6265.csv")
rs4680 = pd.read_csv("/Users/mihriban/Desktop/staj_1/373507-SampleGenotypes-Homo_sapiens_Variation_Sample_rs4680.csv")

rslerlist = [ "rs1815739", "rs4644994", "rs4253778", "rs1042713", "rs8192678", "rs2070744", "rs11549465", "rs1800012", "rs12722", "rs1049434",
        "rs1800795", "rs6265", "rs4680"]

dataframes = [rs1815739, rs4644994, rs4253778, rs1042713, rs8192678, rs2070744, rs11549465, rs1800012, rs12722, rs1049434, rs1800795,
              rs6265, rs4680]

rsler = {str(rs1815739): "rs1815739",str(rs4644994) :"rs4644994", str(rs4253778) :"rs4253778", str(rs1042713) :"rs1042713",
    str(rs8192678) :"rs8192678", str(rs2070744) :"rs2070744", str(rs11549465) :"rs11549465", str(rs1800012) :"rs1800012",
    str(rs12722) :"rs12722", str(rs1049434) :"rs1049434", str(rs1800795) :"rs1800795", str(rs6265) :"rs6265", str(rs4680) :"rs4680"}

def ethnicity(df):
    if os.path.isdir(f"/Users/mihriban/Desktop/staj_1/ethnicity/{rsler[str(df)]}") == False:
        os.mkdir(f"/Users/mihriban/Desktop/staj_1/ethnicity/{rsler[str(df)]}")
    else:
        pass
    #continent = ["AFR","AMR","EAS","EUR","SAS"]
    for i in range(len(df.iloc[:,2])):
        for j in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
            if j in df.iloc[i,2]:
                with open("/Users/mihriban/Desktop/staj_1/ethnicity/{}/{}.csv".format(rsler[str(df)],j),'a+', newline ='') as file:
                    with file:
                        write1 = csv.writer(file)
                        write1.writerow(df.loc[i,:])
for i in dataframes:
    ethnicity(i)









