#Curse of Dimensionality Plot
import matplotlib.pyplot as plt


fig, axs = plt.subplots(2, 3, constrained_layout = True)
fig.set_size_inches(18.5, 10.5)
file = open("datos.txt")

def read_plot1_data():
    title = file.readline()
    vals =[]
    for i in range(11):
        file.read(4)
        vals.append(int(file.readline()))
    file.readline()
    return title, vals

def genPlots():
    x = []
    fig.suptitle("Graficas")
    for i in range(11): x.append(i/10)
    for i in range(2):
        for j in range(3):
            title, y = read_plot1_data()
            axs[i,j].plot(x,y)
            axs[i,j].set_title(title)

genPlots()
for ax in axs.flat:
    ax.set(xlabel='count', ylabel='distance')
    ax.grid(True,linestyle='dotted')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 6500)

plt.savefig('graficas.png')
