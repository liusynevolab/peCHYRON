import pickle
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
from itertools import combinations, product
from scipy.cluster import hierarchy
from scipy.stats import norm
import operator
import numpy as np

LENGTH_CUTOFF = 3 ## inclusive, exclusive
COUNT_CUTOFF = 1.5e-6
RECOVERY_EFFICIENCY = 1

## A class to handle tree construction for comparison purposes
# Tree root is a node of empty character.  Each node holds 
# a dictionary of children, and a value denoting how many sequences
# have the given prefix.
class Tree:
	class Node:	
		def __init__(self, count, char, depth):
			self.char = char
			self.count = count
			self.children = {}
			self.depth = depth
			
	def __init__(self):
		self.root = self.Node(0, '', 0)
		self.length = 0 #max length of a sequence in the tree
		self.nodeList = set()
		
	
	def insertSequence(self, sequence, count):
		curr = self.root #Start at the root
		#curr.count = curr.count + count #Add to the count of this path
		weight = 1/len(sequence)
		curr.count = curr.count + weight
		for char in sequence:
			if char not in curr.children:
				#newNode = self.Node(count, char)
				newNode = self.Node(weight, char, curr.depth + 1)
				curr.children[char] = newNode
				curr = newNode
				self.nodeList.add(newNode)
			else: 
				curr = curr.children[char]
				curr.count = curr.count + weight # count

## A function to get the filenames for all wells.
# Path: The path to the directory containing the folders
# Pattern: The pattern the files take.  Files should be csv headered files
# Returns a dictionary of pandas dataframes of form {name : data}
def getWellFiles(path='', pattern = '*.csv', lengthCutoff = 0, countCutoff = 0):
	# Get list of files to parse in the given path matching the given pattern
	files = glob.glob(os.path.join(path, pattern))
	
	# Parse pattern to allow for name acquisition
	pattern_split = pattern.split("*")
	
	# Make dictionary of data files to return.
	wellData = {}
	for file in files:
		# Strip name of patterns, leading only concatenated wild cards
		name = os.path.basename(file)
		for substring in pattern_split:
			name = name.replace(substring, '')
		
		# Read data in with pandas
		data = pd.read_csv(file, keep_default_na=False)
		countCut = countCutoff * sum(data.iloc[:,2])
		#print(name, countCut)
		data = data[data.iloc[:, 1] > lengthCutoff]

		data = data[data.iloc[:, 2] > countCut]
		wellData[name] = data

	return wellData

## A function to obscure some amount of unique reads from the dataframes.
# wellList: The input of dataFrames
# blindRatio: the proportion of rows to remove
# count: the number of times to try this
# seed: the seed value for the random number generator
# returns a list of wellsLists with obscured data
def obscureData(wellList, blindRatio = 0.75, count = 5, seed=12345):
	obscured = [{} for _ in range(count)]
	rng = np.random.default_rng(seed=seed)
	for name in wellList.keys():
		data = wellList[name]
		for i in range(0, count):
			seed = seed + 10
			rng = np.random.default_rng(seed=seed)
			keep = rng.random(size=len(data))
			newData = data[keep > blindRatio]
			obscured[i][name] = newData
			#print(name, '%7d -> %7d' % (len(data), len(newData)))
	return obscured

## A function to calculate an analogue of generalized jaccard distance between two sets
# This does this by generating trees of inputs, which are then compared and scored
# Based on tree lengths
# x: Tree of sequences
# y: Tree of sequences
# return: a float between 0 to 1 of the jaccard similarity between the two sets
def partialTreeJaccard(xTree, yTree, cutoff=0):
	#Getting the cardinality of the nodes
	xCardinality = 0
	yCardinality = 0
	for node in xTree.nodeList:
		if node.depth > cutoff:
			xCardinality = xCardinality + node.count#/(2**node.depth)
	for node in yTree.nodeList:
		if node.depth > cutoff:
			yCardinality = yCardinality + node.count#/(2**node.depth)
		
	#Comparing trees for partial jaccard distance
	#Explore pairs of nodes in both trees
	exploreList = [(xTree.root, yTree.root)]
	intersectCount = 0
	while len(exploreList)>0:
		#Get node to explore
		currNodePair = exploreList.pop()
		xNode = currNodePair[0]
		yNode = currNodePair[1]
		
		xChildren = xNode.children
		yChildren = yNode.children
		intersect = xChildren.keys() & yChildren.keys()
		for child in intersect:
			exploreList.append((xChildren[child], yChildren[child]))
		if xNode.depth > cutoff:
			intersectCount = intersectCount + min(xNode.count, yNode.count)#/(2**xNode.depth)
	#print(xCardinality, yCardinality, intersectCount)
	return 1 - ((intersectCount)/(xCardinality + yCardinality - intersectCount))

def jaccard(xData, yData):
	xSet = set(xData)
	ySet = set(yData)
	intersect = xSet & ySet
	return 1-len(intersect)/(len(xData) + len(yData) - len(intersect))
	

## Function to generate a dendrogram using Prefix Tree jaccard
# data: A dictionary of dataframes of the data to plot differences of
# cutoff: the depth cutoff for calculating PTJ
# returns: a dataframe consisting of the pairwise distances between entries
def PTJ_dist(data, cutoff=0):
	wells = data.keys()
	wellTrees = {}
	# Generate trees of each well
	for well in wells:
		wellTree = Tree()
		for index, entry in data[well].iterrows():
			wellTree.insertSequence(entry[0], entry[2])
		wellTrees[well] = wellTree

	#Prefix Jaccard Distance
	distDF = pd.DataFrame(columns=[x for x in wells], index=[x for x in wells], dtype='float64')
	for x,y in combinations(wells, 2):
		d = partialTreeJaccard(wellTrees[x], wellTrees[y], cutoff=cutoff)
		#print(x, len(wellDataList[x]), y, len(wellDataList[y]))
		distDF.at[x, y] = d #round(d,1)
		distDF.at[y, x] = d #round(d,1)
		distDF.at[y, y] = 0
		distDF.at[x, x] = 0
	return distDF

## Function to generate a dendrogram using Prefix Tree jaccard
# data: A dictionary of dataframes of the data to plot differences of
# returns: a dataframe consisting of the pairwise distances between entries
def J_dist(data):
	wells = data.keys()
	jacDF = pd.DataFrame(columns=[x for x in wells], index=[x for x in wells], dtype='float64')
	for x,y in combinations(wells, 2):
		d = jaccard(data[x].iloc[:,0], data[y].iloc[:,0])
		jacDF.at[x, y] = d #round(d,1)
		jacDF.at[y, x] = d #round(d,1)
		jacDF.at[y, y] = 0
		jacDF.at[x, x] = 0
	return jacDF


wellDataList = getWellFiles(pattern='*full_signatures_dataframe.csv', lengthCutoff = LENGTH_CUTOFF, countCutoff=COUNT_CUTOFF)
prefixDF = PTJ_dist(wellDataList, 0)
jacDF = J_dist(wellDataList)

#figure, axes = plt.subplots(1, 1)#, figsize=(10, 2))
plt.figure(0)
z = hierarchy.linkage(jacDF, 'average')
hierarchy.set_link_color_palette(['k'])
Index=(['3-1-1','2-1-2','4-1-2','2-2-2','4-2-1','3-1-2','1-1-1','4-2-2','1-2-2','1-2-1','4-1-1','2-2-1','2-1-1'])
dn = hierarchy.dendrogram(z, labels=Index, above_threshold_color='#bbbbbb', orientation='left')
#axes.set_title('Jaccard Dendrogram')
plt.title('Jaccard Dendrogram')
plt.savefig('Jaccard dendrogram')

#figure, axes = plt.subplots(2, 1)#, figsize=(10, 2))
plt.figure(1)
z = hierarchy.linkage(prefixDF, 'average')
hierarchy.set_link_color_palette(['k'])
dn = hierarchy.dendrogram(z, labels=Index, above_threshold_color='#bbbbbb', orientation='left')
#axes.set_title('PTJ Dendrogram')
plt.title('prefix jaccard dendrogram')
plt.savefig('prefix jaccard dendrogram')
#plt.savefig('prefix_jaccard dendrogram')
#plt.show()
