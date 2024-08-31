import pickle
import numpy as np

with open('thickness.pickle', 'rb') as f:
    test = pickle.load(f)
f.close()

print(test[-10:])

with open('thickness_mean.pickle', 'rb') as f:
    test = pickle.load(f)
f.close()

print(test[-10:])


'''
with open('curv_mean.pickle', 'rb') as f:
    test = pickle.load(f)
f.close()

print(test[-10:])

with open('tails_mean.pickle', 'rb') as f:
    test = pickle.load(f)
f.close()

print(test[-10:])

with open('thickness_mean.pickle', 'rb') as f:
    test = pickle.load(f)
f.close()

print(test[-10:])
'''
