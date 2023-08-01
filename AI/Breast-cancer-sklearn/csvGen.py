import pandas as pd

data = pd.read_csv("breast-cancer-wisconsin.data")

data.columns = ["id", "ClumpThick", "UniSize", "UniShape", "MargAd", "SingEpiCelSize", "BareNuc", "BlandChr",
              "NormalNuc", "Mito", "Class"]

data.to_csv("data.csv", index=None, header=True)