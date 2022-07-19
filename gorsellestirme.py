from calisma1part2 import AFR_frequency,AFR_allel_frequency,AMR_frequency, AMR_allel_frequency,EAS_frequency, EAS_allel_frequency,EUR_frequency, EUR_allel_frequency,SAS_frequency, SAS_allel_frequency
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

labels = [ "rs1815739", "rs4644994", "rs4253778", "rs1042713", "rs8192678", "rs2070744", "rs11549465", "rs1800012", "rs12722", "rs1049434","rs1800795", "rs6265", "rs4680"]

unique = ["A", "C", "T", "G"]
uniques = ["A|A","C|C","T|T","G|G","A|C","A|T","A|G","C|T","G|C","T|G"]
#colors = ["cornflowerblue","chocolate","green","gold","purple","darkcyan","saddlebrown","orange","olive"]


def graphing_allel(df_frequency):
    rs_names = []
    rs_values = []
    for i in range(len(df_frequency)):
        rs_names.append(list(df_frequency[i].keys()))
        rs_values.append(list(df_frequency[i].values()))

    handles = []
    def colorfun(j,k):
        color = str()
        if rs_names[j][k] == unique[0]:
            color = "cornflowerblue"
        elif rs_names[j][k] == unique[1]:
            color = "chocolate"
        elif rs_names[j][k] == unique[2]:
            color = "green"
        elif rs_names[j][k] == unique[3]:
            color = "gold"

        return color

    handles.append(mpatches.Patch(color="cornflowerblue", label=unique[0]))
    handles.append(mpatches.Patch(color="chocolate", label=unique[1]))
    handles.append(mpatches.Patch(color="green", label=unique[2]))
    handles.append(mpatches.Patch(color="gold", label=unique[3]))

    #print(handles)
    fig, ax = plt.subplots()

    for j in range(13):

        ax.barh(labels[j], rs_values[j][0], label=rs_names[j][0], color = colorfun(j,0), edgecolor=['none'])
        ax.barh(labels[j], rs_values[j][1], left= rs_values[j][0],label=rs_names[j][1], color= colorfun(j,1), edgecolor=['none'])
        if len(rs_values[j]) == 3:
            ax.barh(labels[j], rs_values[j][2], left=rs_values[j][0]+rs_values[j][1],label=rs_names[j][2], color= colorfun(j,2),edgecolor=['none'])

    #plt.set_ylabel('Scores')
    def suptitle(df_frequency):
        supti = str()
        if df_frequency == AFR_allel_frequency:
            supti = 'AFR allel values'
        elif df_frequency == AMR_allel_frequency:
            supti = 'AMR allel values'
        elif df_frequency == EAS_allel_frequency:
            supti = 'EAS allel values'
        elif df_frequency == EUR_allel_frequency:
            supti = 'EUR allel values'
        elif df_frequency == SAS_allel_frequency:
            supti = 'SAS allel values'
        return supti

    plt.suptitle(suptitle(df_frequency))
    plt.legend(handles = handles, loc="best")

    plt.show()

def graphing(df_frequency):
    rs_names = []
    rs_values = []
    for i in range(len(df_frequency)):
        rs_names.append(list(df_frequency[i].keys()))
        rs_values.append(list(df_frequency[i].values()))

    handles = []
    def colorfun(j,k):
        color = str()
        if rs_names[j][k] == uniques[0]:
            color = "blue"
        elif rs_names[j][k] == uniques[1]:
            color = "chocolate"
        elif rs_names[j][k] == uniques[2]:
            color = "green"
        elif rs_names[j][k] == uniques[3]:
            color = "gold"
        elif rs_names[j][k] == uniques[4]:
            color = "purple"
        elif rs_names[j][k] == uniques[5]:
            color = "darkcyan"
        elif rs_names[j][k] == uniques[6]:
            color = "darkseagreen"
        elif rs_names[j][k] == uniques[7]:
            color = "saddlebrown"
        elif rs_names[j][k] == uniques[8]:
            color = "orange"
        elif rs_names[j][k] == uniques[9]:
            color = "olive"

        return color

    handles.append(mpatches.Patch(color="cornflowerblue", label=uniques[0]))
    handles.append(mpatches.Patch(color="chocolate", label=uniques[1]))
    handles.append(mpatches.Patch(color="green", label=uniques[2]))
    handles.append(mpatches.Patch(color="gold", label=uniques[3]))
    handles.append(mpatches.Patch(color="purple", label=uniques[4]))
    handles.append(mpatches.Patch(color="darkcyan", label=uniques[5]))
    handles.append(mpatches.Patch(color="darkseagreen", label=uniques[6]))
    handles.append(mpatches.Patch(color="saddlebrown", label=uniques[7]))
    handles.append(mpatches.Patch(color="orange", label=uniques[8]))
    handles.append(mpatches.Patch(color="olive", label=uniques[9]))

    #print(handles)
    fig, ax = plt.subplots()

    for j in range(13):

        ax.barh(labels[j], rs_values[j][0], label=rs_names[j][0], color = colorfun(j,0), edgecolor=['none'])
        ax.barh(labels[j], rs_values[j][1], left= rs_values[j][0],label=rs_names[j][1], color= colorfun(j,1), edgecolor=['none'])
        if len(rs_values[j]) == 3:
            ax.barh(labels[j], rs_values[j][2], left=rs_values[j][0]+rs_values[j][1],label=rs_names[j][2], color= colorfun(j,2),edgecolor=['none'])

    #plt.set_ylabel('Scores')
    def suptitle(df_frequency):
        supti = str()
        if df_frequency == AFR_frequency:
            supti = 'AFR values'
        elif df_frequency == AMR_frequency:
            supti = 'AMR values'
        elif df_frequency == EAS_frequency:
            supti = 'EAS values'
        elif df_frequency == EUR_frequency:
            supti = 'EUR values'
        elif df_frequency == SAS_frequency:
            supti = 'SAS values'
        return supti

    plt.suptitle(suptitle(df_frequency))
    plt.legend(handles = handles, loc="best")

    plt.show()

#graphing(AFR_frequency)
#graphing(AMR_frequency)
#graphing(EAS_frequency)   HATA VERİYOR
#graphing(EUR_frequency)
#graphing(SAS_frequency)

#graphing_allel(AFR_allel_frequency)
#graphing_allel(AMR_allel_frequency)
#graphing_allel(EAS_allel_frequency)
graphing_allel(EUR_allel_frequency)
#graphing_allel(SAS_allel_frequency)
