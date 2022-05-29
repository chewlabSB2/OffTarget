#!/usr/bin/env python3

from bitarray import bitarray
import re
from random import choice
from OffTarget.utils import *

class Reference:
	totalCount = 0

	def __init__(self, ID, sequence, adjusted = 0, iteration = 0):
		self.id = ID 
		self.wavelet = Wavelet_Tree(self._replaceUniqueChar(sequence))
		self.indexed = False
		self.process = -1
		self.count = -1
		self.BFbitArray = []
		self.adjusted = adjusted
		self.iteration = iteration

	def _replaceUniqueChar(self, sequence):
		for char in UNIQUE_NT:
			if char in sequence:
				indices = [i.start() for i in re.finditer(char, sequence)]
				for i in indices:
					sequence = sequence[:i] + choice(ALPHABET[char]) + sequence[i + 1:]
		return sequence

	def createBloomFilters(self):
		if self.indexed: return

		def set_bit(value, bit_index):
			return value | (1 << bit_index)

		difference = KMER_LENGTH - 1
		seq = self.wavelet.ReconstructSequence()
		seqLength = self.wavelet.length - difference
		self.BFbitArray = [1 << seqLength for i in range(len(KMER_LIST))]
		for p in range(seqLength):
			mer = seq[p:p+KMER_LENGTH]
			merIndex = KMER_LIST.index(mer)
			value = self.BFbitArray[merIndex]
			self.BFbitArray[merIndex] = set_bit(value, seqLength - p - 1)
		
		self.indexed = True

	def BFdestructor(self):
		del self.BFbitArray	
		self.indexed = False

class Node:
	def __init__(self, value = []):
		self.left  = None #0 
		self.right = None #1
		self.value = bitarray(''.join([str(i) for i in value]))

class Wavelet_Tree():
	
	def __init__(self, sequence):
		sequence = sequence.replace(' ','')
		self.root = None
		self.length = len(sequence)
		self.ALPHABET = list(set(sequence)) 
		self.CharBit = {k:None for k in self.ALPHABET}
		x = {k:sequence.count(k) for k in self.ALPHABET}
		self.nt_count = {k:v for k,v in x.items() if v} 
		self.ConstructNode(sequence, alphabet = self.ALPHABET)
		#self.printTree()

	def ReconstructSequence(self, start = 0, stop = None):
		if not stop: stop = self.length

		assert start >= 0 
		assert stop <= self.length

		sequence = ''
		for i in range(start, stop):
			sequence += self.Access(i)
		return sequence

	def ConstructNode(self, string, node = None, alphabet = None, charbit = '', root = True):
		if len(alphabet) == 1:
			self.CharBit[alphabet[0]] = charbit
			return 

		charbitL = charbitR = charbit
		charbitL += '0'
		charbitR += '1'
		Bitmap = []
		x = {}

		for char in alphabet:
			if char not in list(self.nt_count.keys()): continue
			x[char] = self.nt_count[char]

		alphabet.sort(reverse = True)
		keys = alphabet
		mid = len(keys)//2
		SplitCharL = keys[:mid]
		SplitCharR = keys[mid:]
		Sleft = Sright = ''

		for c in list(string):
			if c in SplitCharL:
				Sleft += c
				Bitmap.append(0)
			else:
				Sright += c
				Bitmap.append(1)

		newNode = Node(Bitmap)
		if not self.root:
			self.root = newNode

		if charbit: 
			current_charbit = charbit[-1]
			if current_charbit == '0':
				node.left = newNode
			elif current_charbit == '1':
				node.right = newNode

		self.ConstructNode(Sleft, newNode, SplitCharL, charbitL, root = False)
		self.ConstructNode(Sright, newNode, SplitCharR, charbitR, root = False)
		if not root: return
		
		self.CharBit = {k:v for k,v in self.CharBit.items() if v}
		self.CharBit_inv = {v:k for k, v in self.CharBit.items()}

	def printTree(self):
		if self.root is not None:
			self._printTree(self.root)

	def _printTree(self, node, left = True):
		if node is not None:
			self._printTree(node.left, left)
			print (left, node.value.__sizeof__(), str(node.value) + ' ')
			self._printTree(node.right, False)

	def _BinaryRank(self, c, p, node):
		return node.value[:p].count(c)

	def Rank(self, c, p):
		'''	
		c = Character
		p = position/ offset
		'''
		if p < 0:
			return 0
		charbit = list(self.CharBit[c])
		max_counter = len(charbit)
		i = 0
		node = None
		while True: 
			if not node:
				node = self.root
			else:
				if current_charbit == 0:
					node = node.left
				elif current_charbit == 1:
					node = node.right
			
			current_charbit = int(charbit[i])
			p = self._BinaryRank(current_charbit, p, node)
			i += 1
			if i == len(charbit):
				break
				
		return p

	def _BinaryAccess(self, node, o):
		return node.value[o]

	def Access(self, offset):
		charbit = ''
		node = None
		while True:
			if not node:
				node = self.root
			o = int(self._BinaryAccess(node, offset))
			charbit += str(o)

			if charbit in list(self.CharBit_inv.keys()):
				return self.CharBit_inv[charbit]

			r = int(self._BinaryRank(o, offset, node))
			offset = r 
			if o == 0:
				node = node.left
			elif o == 1:
				node = node.right
	
	def Select(self):
		pass