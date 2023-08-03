import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC

import pickle

# Read the data, add the column names and save it to a csv.     
data = pd.read_csv("breast-cancer-wisconsin.data")

data.columns = ["id", "ClumpThick", "UniSize", "UniShape", "MargAd", "SingEpiCelSize", "BareNuc", "BlandChr",
              "NormalNuc", "Mito", "Class"]

data.to_csv("data.csv", index=None, header=True)

# Read the data and preprocess it.
data = pd.read_csv("data.csv")

data.drop(['id'], inplace=True, axis = 1)
data.replace('?', -99999, inplace = True)

data["Class"] = data["Class"].map(lambda x: 1 if x == 4 else 0)

print(data.head())

# Select x and y data for the models
X = np.array(data.drop(["Class"], axis = 1))
y = np.array(data["Class"])

# Training and testing the models

[X_train, X_test, y_train, y_test] = train_test_split(X, y, test_size = 0.1, random_state = 0)

# # SVC classifier
Classifier = SVC(kernel = "linear")
model = Classifier.fit(X_train, y_train)
accu = model.score(X_test, y_test)
print("Accuracy of SVC: ", accu)

# # Logistic Regression
Classifier2 = LogisticRegression(solver = 'liblinear')
model2 = Classifier2.fit(X_train, y_train)
accu2 = model2.score(X_test, y_test)
print("Accuracy of Logistic Regression: ", accu2)

# Save the model using pickle
pickle.dump(model, open("LogisticRegressionBreastCancerSKL.m", "wb"))

loaded_model = pickle.load(open("LogisticRegressionBreastCancerSKL.m", "rb"))
accu = loaded_model.score(X_test, y_test)
print("Accuracy of loaded model: ", accu)

# Predict the tumour based on the characteristics
classes = ["Benign", "Malignant"]
sample = np.array([[8,10,10,8,7,10,9,7,1]]) # Characteristics based on the data, 
                                            # change these numbers for different predictions
result = loaded_model.predict(sample)
print(classes[int(result)])