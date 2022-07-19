from calisma1part2 import AFR, AMR, EAS, EUR, SAS,heterovals


Power = {"rs1815739": {'T|T':0,'C|T':1,'C|C':2},
    "rs4253778": {'G|G':0,'G|C':1,'C|C':2},
    "rs1042713": {'A|A':0,'A|G':1,'G|G':2},
    "rs8192678": {'C|C':0,'C|T':1,'T|T':2},
    "rs11549465": {'C|C':0,'C|T':1,'T|T':2} }

Endurance = {"rs1815739" : {'T|T':2,'C|T':1,'C|C':0},
    "rs4253778": {'G|G':0,'G|C':1,'C|C':2},
    "rs1042713": {'A|A':2,'A|G':1,'G|G':0},
    "rs8192678": {'C|C':2,'C|T':1,'T|T':0},
    "rs11549465": {'C|C':2,'C|T':1,'T|T':0}}

Injury_Risk = {"rs1800012" : {'C|C':2,'A|C':1,'A|A':0},
    "rs12722": {'C|C':2,'C|T':1,'T|T':0},
    "rs1049434": {'A|A':2,'A|T':1,'T|T':0},
    "rs1800795": {'G|G':0,'G|C':1,'C|C':2} }

Oxygen_Capacity={"rs11549465": {'C|C':2,'C|T':1,'T|T':0} }

Motor_Skills = { "rs6265": {'C|C':2,'C|T':1,'T|T':0},
    "rs4680": {'A|A':2,'A|G':1,'G|G':0} }


rslerdict = { "rs1815739" :0, "rs4644994":1, "rs4253778":2, "rs1042713":3, "rs8192678":4, "rs2070744":5, "rs11549465":6,
              "rs1800012":7,"rs12722":8, "rs1049434":9,"rs1800795":10, "rs6265":11, "rs4680":12}

G_i = 0
n = 0

def df_change(df):

    for i in range(len(df)):
        for j in range(13):
            for a,b in heterovals:
                if a == df.iloc[i,j] or b == df.iloc[i,j]:
                    df.iloc[i,j] = a
    return df


AFR = df_change(AFR)
AMR = df_change(AMR)
EAS = df_change(EAS)
EUR = df_change(EUR)
SAS = df_change(SAS)


def gen_score(df):

    power_scores = []
    endurance_scores = []
    injury_scores = []
    oxygen_scores = []
    motor_scores = []

    power_list = [0,2,3,4,6]
    endurance_list = [0,2,3,4,6]
    injury_list = [7,8,9,10]
    oxygen_list = [6]
    motor_list = [11,12]

    df_power = df.iloc[:, power_list]
    df_power = df_power.replace(Power)
    power_scores.append((df_power.sum(axis=1)/(2 * len(power_list)))*100)

    df_endurance = df.iloc[:, endurance_list]
    df_endurance = df_endurance.replace(Endurance)
    endurance_scores.append((df_endurance.sum(axis=1) / (2 * len(endurance_list))) * 100)

    df_injury = df.iloc[:, injury_list]
    df_injury = df_injury.replace(Injury_Risk)
    injury_scores.append((df_injury.sum(axis=1) / (2 * len(injury_list))) * 100)

    df_oxygen = df.iloc[:, oxygen_list]
    df_oxygen = df_oxygen.replace(Oxygen_Capacity)
    oxygen_scores.append((df_oxygen.sum(axis=1) / (2 * len(oxygen_list))) * 100)

    df_motor = df.iloc[:, motor_list]
    df_motor = df_motor.replace(Motor_Skills)
    motor_scores.append((df_motor.sum(axis=1) / (2 * len(motor_list))) * 100)

    return power_scores, endurance_scores, injury_scores, oxygen_scores, motor_scores


power_scores_AFR, endurance_scores_AFR, injury_scores_AFR, oxygen_scores_AFR, motor_scores_AFR = gen_score(AFR)
power_scores_AMR, endurance_scores_AMR, injury_scores_AMR, oxygen_scores_AMR, motor_scores_AMR = gen_score(AMR)
power_scores_EAS, endurance_scores_EAS, injury_scores_EAS, oxygen_scores_EAS, motor_scores_EAS = gen_score(EAS)
power_scores_EUR, endurance_scores_EUR, injury_scores_EUR, oxygen_scores_EUR, motor_scores_EUR = gen_score(EUR)
power_scores_SAS, endurance_scores_SAS, injury_scores_SAS, oxygen_scores_SAS, motor_scores_SAS = gen_score(SAS)




