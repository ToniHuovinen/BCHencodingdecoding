# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 20:03:01 2018

@author: Toni Huovinen
"""


import numpy as np

# Generator Polynomial and word for encoding in binary form
Gbin = np.array([1,0,0,0,1,0,1,1,1])
mbin = np.array([1,1,1,0,0,0,1])


# Encoding
def binary_multiplication(G, m):
    fullbase = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    tempbase = np.array([0,0,0,0,0,0,0,0,0])

    for i in range(len(m[:])):
        for j in range(len(G[:])):
            temp = tempbase[:]
            temp[8-j] = m[6-i] * G[8-j]
        
            fullbase[(14-i)-j] = temp[8-j] + fullbase[(14-i)-j]


    fullbase = np.mod(fullbase, 2)
    
    return fullbase

C = binary_multiplication(Gbin, mbin)

print("Word to be encoded:")
print(mbin)
print("Encoded it is:")
print(C)


#===================================================================


# Decoding
# Received word (includes two errors)
r = np.array([1,1,1,0,1,0,0,1,0,0,1,1,1,1,0])

field = np.array([[1,0,0,0], # 1
                  [0,1,0,0], # a
                  [0,0,1,0], # a**2
                  [0,0,0,1], # a**3
                  [1,1,0,0], # 1 + a
                  [0,1,1,0], # a + a**2
                  [0,0,1,1], # a**2 + a**3
                  [1,1,0,1], # 1 + a + a**3
                  [1,0,1,0], # 1 + a**2
                  [0,1,0,1], # a + a**3
                  [1,1,1,0], # 1 + a + a**2
                  [0,1,1,1], # a + a**2 + a**3
                  [1,1,1,1], # 1 + a + a**2 + a**3
                  [1,0,1,1], # 1 + a**2 + a**3
                  [1,0,0,1]]) # 1 + a**3

betachart = np.array([[5,10], # Match for 0,1 is not on the list since we want to fix
                      [6,13], # two bits
                      [11,12],
                      [99, 99],
                      [7,9],
                      [2,8],
                      [99, 99],
                      [99, 99],
                      [3, 14],
                      [99, 99],
                      [1, 4],
                      [99, 99],
                      [99, 99],
                      [99, 99],
                      [99, 99]])
                         
    
# Syndromes
# Base
S1 = np.array([0,0,0,0])
S3 = np.array([0,0,0,0])

# Polynomial degrees within the received word
degrees = np.array([])
for i in range(len(r[:])):
    if r[i] == 1:
        degrees = np.append(degrees, i)


# Function for calculating the syndromes      
def syndrome_calculation(degrees, fields, syndrome):
    
    for i in range(len(degrees[:])):
        for j in range(15):
            if degrees[i] == j:
                syndrome = np.add(syndrome, fields[j])

    return syndrome


# Syndrome S1     
S1 = syndrome_calculation(degrees, field, S1) 
S1 = np.mod(S1, 2)


# Syndrome S3
for i in range(len(degrees[:])):
    degrees[i] = degrees[i] * 3
    
degrees = np.mod(degrees, 15)

S3 = syndrome_calculation(degrees, field, S3)
S3 = np.mod(S3, 2)


# Degrees of syndromes
def syndrome_degree(field, syndrome):
    for i in range(len(field[:])):
        if np.array_equal(syndrome[:], field[i]):
            degree = i;
    return degree


degreeS1 = syndrome_degree(field, S1)
degreeS3 = syndrome_degree(field, S3)

# Find the match based on beta chart
def beta_match(S1_degree, S3_degree, match_chart):
    S1_3 = S1_degree * 3
    beta = S3_degree - S1_3
    
    if beta < 0:
        beta = beta + 15
        
    matches = match_chart[beta]
    
    for i in range(2):
        matches[i] = matches[i] + S1_degree
        
    matches = np.mod(matches, 15)
    
    return matches


matches = beta_match(degreeS1, degreeS3, betachart)

def print_decoded(Gbin, received_word, match):
    errorbase = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    errorbase[match[0]] = 1
    errorbase[match[1]] = 1


    Cstar = np.add(received_word, errorbase)
    Cstar = np.mod(Cstar, 2)

    mstar = np.polydiv(Cstar, Gbin)
    mstar = np.mod(mstar, 2)
    
    print("Decoded word is:")
    print(mstar[0])

# Result
print_decoded(Gbin, r, matches)



