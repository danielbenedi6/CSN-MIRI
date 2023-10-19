import pandas
import matplotlib.pyplot as plt
import seaborn
languages = ['Arabic', 'Basque', 'Catalan', 'Chinese', 'Czech', 'English', 'Greek', 'Hungarian', 'Italian', 'Turkish']
types = ["binomial", "switching"]

df = [
    ("Arabic", "Graph", 0.32597),
    ("Basque", "Graph", 0.26981),
    ("Catalan", "Graph", 0.34222),
    ("Chinese", "Graph", 0.32564),
    ("Czech", "Graph", 0.30661),
    ("English", "Graph", 0.34381),
    ("Greek", "Graph", 0.31403),
    ("Hungarian", "Graph", 0.28971),
    ("Italian", "Graph", 0.32763),
    ("Turkish", "Graph", 0.36075)
]

for lan in languages:
    for type in ["binomial", "switching"]:
            f = open("results/%s_%s.csv"%(lan,type),"r")
            lines = f.readlines()
            df += [(lan, type, float(line)) for line in lines]

df = pandas.DataFrame(df, columns=["Language", "Model", "Closeness"])

seaborn.stripplot(data=df[df["Model"] != "Graph"],x="Closeness",y="Language", hue="Model")
seaborn.boxplot(
                data=df[df["Model"] == "Graph"],
                x="Closeness",
                y="Language",
                showmeans=True,
                meanline=True,
                meanprops={'color': 'k', 'ls': '-', 'lw': 2},
                medianprops={'visible': False},
                whiskerprops={'visible': False},
                zorder=10
                )
plt.title("Comparison of closeness centrality of the models and the input graph")
plt.savefig("./closeness_distribution.pdf")
plt.show()
plt.cla()
plt.clf()

validation = pandas.read_csv("./validation.csv", sep=";")

seaborn.stripplot(data=validation[validation["Type"] != "Exact"],x="Closeness",y="Model")
seaborn.boxplot(
                data=validation[validation["Type"] == "Exact"],
                x="Closeness",
                y="Model",
                showmeans=True,
                meanline=True,
                meanprops={'color': 'k', 'ls': '-', 'lw': 2},
                medianprops={'visible': False},
                whiskerprops={'visible': False},
                zorder=10
                )
plt.title("Comparison of montecarlo approximation with the exact evaluation")
plt.savefig("./closeness_validation.pdf")
plt.show()
plt.cla()
plt.clf()

