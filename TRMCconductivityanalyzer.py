#imports
import pandas as pd
import matplotlib.pyplot as plt
import tkinter as tk
import statsmodels.api as sm
from tkinter import filedialog

#gets data from CSV with conductance data
root = tk.Tk()
root.withdraw()
FN = filedialog.askopenfilename()
data1 = pd.read_csv(FN,skiprows = 13, usecols = [0,2])
data1.columns= ['Time(ns)','Conductance(S)']
data1['Time(ns)'] = data1['Time(ns)'].apply(lambda x: x*1e9)

#Get predictions
X = data1['Time(ns)'] ## X usually means our input variables (or independent variables)
y = data1['Conductance(S)'] ## Y usually means our output/dependent variable
X = sm.add_constant(X) ## adds an intercept (beta_0) to our model
model = sm.OLS(y, X).fit() ## sm.OLS(output, input)
predictions = model.predict(X)

#Add predictions to dataframe
data1['Linear Regression']= predictions

#plot
ax = plt.gca()
data1.plot(kind='scatter',logy=True, x='Time(ns)',y='Conductance(S)',ax=ax, s=.5)
data1.plot(kind='scatter',logy=True, x='Time(ns)',y='Linear Regression',ax=ax,color='red',s=.5)
ax.set_ylabel("Conductance(S)")
plt.autoscale(enable=True, axis='y')

#finds minimum and index
Min = data1['Conductance(S)'].min()
val_mask = data1['Conductance(S)'] == Min
ind = data1['Conductance(S)'][val_mask].index

# finds value of lin regression at index
pt = data1['Linear Regression'][ind].values

#Finds and displays difference
Diff=pt-Min
print('The difference between the minimum and linear regression is',Diff[0],'volts for', FN)

#Display dataframe. Just uncomment if you want it.
data1.head()