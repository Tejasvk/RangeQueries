# Error as a function of b and CI.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, font_manager
import matplotlib as mpl
import sys

mpl.rc('text', usetex=True)

mpl.rcParams['axes.linewidth'] = 3  # set the value globally
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']


d = 8
olh = False
#df = pd.read_csv("../vary_b_ci_"+str(d)+".csv")
#df = pd.read_csv("../vary_b_ci_"+str(d)+"_prefix.csv")
cols = ["length","base","OUE","TreeOUE","TreeHRR","TreeOUECI","TreeHRRCI","HaarHRR"]

#df = pd.read_csv("../vary_b_ci_"+str(d)+".csv",usecols=cols)
#df = pd.read_csv("../df_22.csv",usecols=cols)
df = pd.read_csv("../vary_b_ci_"+str(d)+".csv")
#df = pd.read_csv("../vary_b_ci_"+str(d)+"_prefix.csv")

df = df[1:]

df.rename(columns={"InputRR":"OUE","TreeRR":"TreeOUE","TreeHT": "TreeHRR", "TreeHTCI": "TreeHRRCI",
                   "TreeRRCI": "TreeOUECI"}, inplace=True)
#population = "$2^{"+str(int(np.log2(df.iloc[5]["population"])))+"}$"
population = "$2^{"+str(26)+"}$"

print (population)
#df.drop(["population", "count"], inplace=True, axis=1)
print(df.columns)


cols = df.columns.tolist()
lengths = df["length"].unique()
#print (lengths)
i = 0

if d == 8:
    limindex = {0: (0.0, 0.004), 1: (0.0, 0.007), 2: (0.0, 0.007), 3: (0.0, 0.007),4: (0.0, 0.007)}
    perc_array = [1,15,100,200, np.power(2,d)-2]#np.arange(0,100000,15000)
elif d == 16:
    perc_array = [1,20,1000,10000, np.power(2,d)-2]#np.arange(0,100000,15000)
    limindex = {0: (0.0, 0.004), 1: (0.0, 0.007), 2: (0.0, 0.007), 3: (0.0, 0.007),4: (0.0, 0.02)}
elif d == 18:
    perc_array = [1,32 ,1000, 10000, np.power(2,d)-2]  # np.arange(0,100000,15000)
    limindex = {0: (0.0, 0.004), 1: (0.0, 0.009), 2: (0.0, 0.009), 3: (0.0, 0.009),4: (0.0, 0.04)}
elif d == 20: 
    perc_array = [1, 42, 1000, 100000, np.power(2,d)-2]  # np.arange(0,100000,15000)
    limindex = {0: (0.0, 0.004), 1: (0.0, 0.009), 2: (0.0, 0.009), 3: (0.0, 0.009),4: (0.0, 0.09)}
elif d == 22: 
    perc_array = [1, 42, 1000, 100000, np.power(2,d)-2]  # np.arange(0,100000,15000)
    limindex = {0: (0.0, 0.004), 1: (0.0, 0.012), 2: (0.0, 0.012), 3: (0.0, 0.012),4: (0.0, 0.12)}
    


df = df[df['length'].isin(perc_array)]

dfmean = pd.DataFrame(columns=cols)
dfstd = pd.DataFrame(columns=cols)
dfgrp = df.groupby(["length", "base"])
i = 0
for g, r in dfgrp:

    dfmean.loc[i] = r.mean()
    dfstd.loc[i] = r.std()
    i += 1

dfstd["length"] = dfmean["length"]
dfstd["base"] = dfmean["base"]
# print(dfmean)


dfmeangrp = dfmean.groupby(["length"])
i = 0
dfarray = []
for g, r in dfmeangrp:
    #print ("----------------")
    #print (r)
    dfarray.append([r, g])
    i += 1
dfstdgrp = dfstd.groupby(["length"])
#print (dfstd)
i = 0

#print (dfstd)
for g, r in dfstdgrp:
    #print ("================")
    #print (r)
    dfarray[i].append(r)
    #print (dfarray[i][2].round(5))
    i += 1
# sys.exit(0)
dffinal = pd.DataFrame(columns=[cols], dtype="float64")
dfstd = pd.DataFrame(columns=[cols], dtype="float64")
rows = 1
columns = len(perc_array)
#print (len(dfarray))
cnt = 0
array = np.array(range(0, len(dfarray))).reshape((columns, rows))
fig, axes = plt.subplots(nrows=rows, ncols=columns,
                         figsize=(17, 3), sharey=False)

sizeOfFontX = 18
ticks_fontX = font_manager.FontProperties(
    style='normal',  size=sizeOfFontX, weight='bold', stretch='normal')
sizeOfFontY = 16
ticks_fontY = font_manager.FontProperties(
    style='normal',  size=sizeOfFontY, weight='bold', stretch='normal')

#print (array)
cnt = 0
#flat = "HRR"
flat = "OUE"

if olh == True:
    Allmechanisms = [flat, "TreeOUE", "TreeOUECI", "TreeHRR", "TreeHRRCI", "HaarHRR" ,"TreeOLH","TreeOLHCI"]
    clrs = ["red", "maroon", "maroon", "darkorange", "darkorange", "green","black","black"]  # ,"slateblue"]
else:
    Allmechanisms = [flat, "TreeOUE", "TreeOUECI", "TreeHRR", "TreeHRRCI", "HaarHRR" ]#,"TreeOLH","TreeOLHCI"]
    clrs = ["red", "maroon", "maroon", "darkorange", "darkorange", "green"]#,"black","black"]

#limindex = {0:(0.0,0.004),1:(0.0,0.01),2:(0.0,0.01),3:(0.0,0.01)}

for i in range(len(array)):
    try:
        #print ("---------------------------")
        #print (int(dfarray[cnt][1]))
        axes[i].set_title("N=" + str(population) + ", r="+str(int(dfarray[cnt][1])) +   ", $D=2^{"+str(d)+"}$", fontweight="bold", size=14)
        dfbase = dfarray[cnt][0]  # .copy()
        #print (dfbase[Allmechanisms])
        indlist = dfbase.base.astype("int").tolist()
        #print (indlist)
        dfbase.index = map(lambda x: "$2^{"+str(int(np.log2(x)))+"}$", indlist)

        dfbase2 = dfbase[Allmechanisms].copy()
        dfbasestd = dfarray[cnt][2]  # .copy()
        dfbasestd.index = dfbase.index
        dfbasestd2 = dfbasestd[Allmechanisms].copy()
        #print (dfarray[cnt][1])
        ind = "$2^{"+str(d)+"}$"
        if olh == True:
            dfbase2.loc[ind] = np.zeros(8)
            dfbasestd2.loc[ind] = np.zeros(8)
        else:
            dfbase2.loc[ind] = np.zeros(6)
            dfbasestd2.loc[ind] = np.zeros(6)
        dfbase2[flat] = df[df["length"] == perc_array[i]][flat].mean()
        dfbasestd2[flat] = df[df["length"] == perc_array[i]][flat].std()
        dfbasestd2[flat][0:-1] = 0.0
        dfbase2[flat][0:-1] = 0.0

        dfbase2["HaarHRR"] = df[df["length"]  == perc_array[i]]["HaarHRR"].mean()
        dfbasestd2["HaarHRR"] = df[df["length"] == perc_array[i]]["HaarHRR"].std()
        dfbasestd2['HaarHRR'][1:] = 0.0
        dfbase2['HaarHRR'][1:] = 0.0

        dfbase2.plot.bar(ax=axes[i], color=clrs, rot=0,
                         yerr=dfbasestd2, ylim=limindex[i], width=0.8,legend=False)
        if (cnt == 0):
            axes[i].legend(loc="best", prop={'size': 11, "weight": "bold"})
        #axes[i].set_xlabel("B",fontweight="bold",size=20)

        for label in axes[i].get_xticklabels():
            label.set_fontproperties(ticks_fontX)
        for label in axes[i].get_yticklabels():
            label.set_fontproperties(ticks_fontY)

        cnt += 1
    except:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)


plt.subplots_adjust(wspace=0.4, hspace=0.6, left=0.1,
                    top=0.9, right=0.95, bottom=0.2)
plt.show()

#plt.savefig("d_"+str(d)+"_ci")
plt.close()
