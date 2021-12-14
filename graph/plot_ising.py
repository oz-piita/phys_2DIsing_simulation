import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path = "."
df_m = pd.read_csv(path+'/result/metropolis.csv')
df_g = pd.read_csv(path+'/result/gibbs.csv')
df_a = pd.read_csv(path+'/result/analitic.csv')

# Energy
x = df_m["逆温度"].tolist()
y = df_m['内部エネルギー'].tolist()
err = df_m['内部エネルギー標準誤差'].tolist()
plt.errorbar(x,y,yerr=err,label="metropolis")

plt.xlabel('beta')
plt.ylabel('energy')
plt.legend()
plt.savefig(path+'/graph/met_ene.png',dpi=144)

y = df_g['内部エネルギー'].tolist()
err = df_g['内部エネルギー標準誤差'].tolist()
plt.errorbar(x,y,yerr=err,label="gibbs")

plt.xlabel('beta')
plt.ylabel('energy')
plt.legend()
plt.savefig(path+'/graph/gib_ene.jpg',dpi=144)

y = df_a['内部エネルギー'].tolist()
plt.plot(x,y,marker='o',label="analitical solution")

plt.xlabel('beta')
plt.ylabel('energy')
plt.legend()
plt.savefig(path+'/graph/ana_ene.jpg',dpi=144)

# Magnetization
x = df_m["逆温度"].tolist()
y = df_m['自発磁化'].tolist()
err = df_m['磁化標準誤差'].tolist()
plt.errorbar(x,y,yerr=err,label="metropolis")

plt.xlabel('beta')
plt.ylabel('magnetization')
plt.legend()
plt.savefig(path+'/graph/met_mag.jpg',dpi=144)

y = df_g['自発磁化'].tolist()
err = df_g['磁化標準誤差'].tolist()
plt.errorbar(x,y,yerr=err,label="gibbs")

plt.xlabel('beta')
plt.ylabel('magnetization')
plt.legend()
plt.savefig(path+'/graph/gib_mag.jpg',dpi=144)

y = df_a['自発磁化'].tolist()
plt.plot(x,y,marker='o',label="analitical solution")

plt.xlabel('beta')
plt.ylabel('magnetization')
plt.legend()
plt.savefig(path+'/graph/ana_mag.jpg',dpi=144)

# Capacity
x = df_m["逆温度"].tolist()
y = df_m['比熱'].tolist()
err = df_m['比熱標準誤差'].tolist()
plt.errorbar(x,y,yerr=err,label="metropolis")

plt.xlabel('beta')
plt.ylabel('Capacity')
plt.legend()
plt.savefig(path+'/graph/met_cap.png',dpi=144)

y = df_g['比熱'].tolist()
err = df_g['比熱標準誤差'].tolist()
plt.errorbar(x,y,yerr=err,label="gibbs")

plt.xlabel('beta')
plt.ylabel('Capacity')
plt.legend()
plt.savefig(path+'/graph/gib_cap.jpg',dpi=144)

y = df_a['比熱'].tolist()
plt.plot(x,y,marker='o',label="analitical solution")

plt.xlabel('beta')
plt.ylabel('Capacity')
plt.legend()
plt.savefig(path+'/graph/ana_cap.jpg',dpi=144)