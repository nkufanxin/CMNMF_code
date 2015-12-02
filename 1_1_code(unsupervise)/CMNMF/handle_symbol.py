filename = 'symbol_idx.txt'
filename2 = 'entrez.txt'
   
f2 = open(filename2)         

symbol_set = []         
import fileinput 
for line in fileinput.input(filename):
    symbol_set.append(int(line)-1)
symbols = []
for line in fileinput.input(filename2):
    symbols.append((line))
# File: readline-example-1.py
 
symboles_1274 = []
for item in symbol_set:
	symboles_1274.append(symbols[item]) 
#print symboles_1274

file_object = open('symbols_1274.txt', 'w')
file_object.writelines(symboles_1274)