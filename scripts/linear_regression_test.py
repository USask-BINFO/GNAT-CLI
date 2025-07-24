#! /usr/bin/env python3

from sklearn.linear_model import LinearRegression
from scipy.stats import norm
from math import sqrt
import pandas as pd
import json

data = {}
with open("summary(4).json", "r") as f:
    data = json.load(f)

for query in data.get('query_names'):
    x1, x2, y, path = [], [], [], []

    for plot_item in data['query_names'][query].get('plot_paths'):
        counts = plot_item['counts']
        if counts:
            total, labelled, gns = counts[0], counts[1], counts[2]

            x1.append(total)
            x2.append(gns)
            y.append(labelled)
            path.append(plot_item['path'].split('/')[-1])

df = pd.DataFrame({'path': path, 'x1': x1, 'x2': x2, 'y': y})

model = LinearRegression()

X, y = df[['x1', 'x2']], df.y

model.fit(X, y)

print(model.intercept_, model.coef_, model.score(X, y))

pred = model.predict(X)
res = y - pred

df['residuals'] = res

# Get the Z score of the residuals
# X - mean / stdev
# The mean should be ~0 since we're dealing with residuals (should be an ~equal number of values above and below the linear regression line)
# Check
N = len(res)
t = sum(res)
print(t/N)

# Variance is (1/N) * sum((X-u)^2) since u is ~0, we can just do (1/N) * sum(X^2)
# Std Dev is sqrt(var)
squared_res = res ** 2
print(res)
variance = (1/len(res)) * sum(squared_res)
print(squared_res)
print(variance)

stdev = sqrt(variance)
print(stdev)

# Since mean approx == 0, Z-score can just be X / stdev
Z = res / stdev

print(Z)

print(df[Z > 1])

df['z_scores'] = Z
df['p_values'] = 1 - norm.cdf(Z)

for i,path in enumerate(df.path):
    print(i, path)

print()

print(df[df['p_values'] <= 0.15].sort_values('p_values'))
# for i in df[Z > 1]:
#     print(i)

print(len(df[res > 0].path))
print(len(y))


high_labels = df[df['p_values'] <= 0.15].sort_values('p_values')

print()
print()

print(high_labels[['path', 'p_values']])

for row in high_labels.itertuples(index=False):
    print(row.path, f"{row.p_values:.2f}", round(row.p_values, 2))