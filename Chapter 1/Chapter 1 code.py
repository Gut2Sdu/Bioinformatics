# In[0] 生物信息学算法

import numpy as np
import matplotlib.pyplot as plt
import time
from collections import defaultdict


# with open("C:/Users/Lenovo/Desktop/大肠杆菌基因序列.txt", "r") as f:
#     data = f.read() 
# data.replace(" ","")
# lines = open('C:/Users/Lenovo/Desktop/大肠杆菌基因序列.txt').readlines()
# data= lines[1:len(lines)]
# data = "".join(data)
# data=data.replace('\n', '')


with open("C:/Users/Lenovo/Desktop/E-coli.txt", "r") as f:
    data = f.read()


# In[1] DNA复制从基因组的何处开始

def PatternCount(Text, Pattern): 
    count = 0 
    for i in range(0, len(Text) - len(Pattern) + 1):
        if text[i : len(Pattern) + i] == Pattern:
            count = count + 1 
    return count

# text = "ACAACTATGCATACTATCGGGAACTATCCT"
# pattern = "ACTAT"
# number=PatternCount(text, pattern)


def FrequentWords(Text, k):
    FrequentPatterns = []
    count=[]
    for i in range(0, (len(Text) - k + 1)):
        Pattern=Text[i : i+k] 
        count.append(PatternCount(Text, Pattern))
        maxCount=max(count)
    for i in range(0, (len(Text) - k + 1)):
        if count[i]==maxCount:
            FrequentPatterns.append(Text[i : i+k])
    FrequentPatterns=list(set(FrequentPatterns))
    return FrequentPatterns
    # frequent_patterns = []
    # count = {}
    # for i in range(len(Text) - k):
    #     pattern = Text[i : i + k]
    #     count[i] = PatternCount(Text, pattern)
    #     max_count = max(count.values())
    # for position in range(len(Text) - k):
    #     if count[position] == max_count:
    #         frequent_patterns.append(Text[position : position + k])
    # return frequent_patterns

# text = "ACTGACTCCCACCCC"
#tic = time.time()
# pattern=FrequentWords(text, 3)
#toc = time.time()
#shijian = toc-tic
# print(shijian)

def ReverseComplement(Text):
    complement = []
    for i in Text:
        if i == "A":
            complement.append("T")
        elif i == "T":
            complement.append("A")
        elif i == "G":
            complement.append("C")
        elif i == "C":
            complement.append("G")
        else:
            print("Not ACTG!")
    reverse_complement = list(reversed(complement))
    return "".join(reverse_complement)

# pattern=ReverseComplement("AAAACCCGGT")


def HammingDistance(p, q):
    d = 0
    for p, q in zip(p, q):
        if p!= q:
            d += 1
    return d

#text1="ABCD"
#text2="AACB"
#number=HammingDistance(text1,text2)


def ApproximatePatternCount(Text,Pattern,d):
    count=0
    k=len(Pattern)
    for i in range(0, len(Text) - len(Pattern) + 1):
        Pattern_temp=Text[i : i+k]
        if HammingDistance(Pattern,Pattern_temp)<=d:
            count=count+1
    return count

# text="AACAAGCATAAACATTAAAGAG"
# pattern="AAAAA"
# d=1
# number=ApproximatePatternCount(text,pattern,d)


def Prefix(Pattern):
    return Pattern[:-1] #去掉尾部一个字符

# #pattern="ACGT"
# #result=Prefix(pattern)

def LastSymbol(Pattern):
    return Pattern[-1]

# #pattern="ACGT"
# #result=LastSymbol(pattern)

def SymbolToNumber(Symbol):
    if Symbol=='A':
        return 0
    elif Symbol=='C':
        return 1
    elif Symbol=='G':
        return 2
    elif Symbol=='T':
          return 3

# #symbol='C'
# #number=SymbolToNumber(symbol)

def PatternToNumber(Pattern):
    if Pattern.strip()=='':
        return 0
    else:
        symbol=LastSymbol(Pattern)
        prefix=Prefix(Pattern)
        return 4*PatternToNumber(prefix)+SymbolToNumber(symbol)
    # Pattern=Pattern[::-1]#将dna片段反转，利于结果计算
    # result=0
    # i=0
    # while i<len(Pattern):
    #     if Pattern[i]=='A':
    #         pass
    #     elif Pattern[i]=='C':
    #         result+=4**i
    #     elif Pattern[i]=='G':
    #         result+=2*4**i
    #     else:
    #         result+=3*4**i
    #     i+=1
    # return result

# #pattern="CG"
# #number=PatternToNumber(pattern)
    
def NumberToPattern(i,k):
    Pattern=[]
    while k>0:
        if i//4**(k-1)==3:
            Pattern.append("T")
        elif i//4**(k-1)==2:
            Pattern.append("G")
        elif i//4**(k-1)==1:
            Pattern.append("C")
        elif i//4**(k-1)==0:
            Pattern.append("A")
        i=i%4**(k-1)
        k=k-1
    Pattern= ''.join(Pattern)
    return Pattern

# #pattern=NumberToPattern(6,2)


def ComputingFrequencies(Text,k):
    FrequencyArray=[]
    for i in range(0, np.power(4,k)):
        FrequencyArray.append(0)
    for i in range(0, (len(Text) - k + 1)):
        Pattern_temp=Text[i : i+k]
        j=PatternToNumber(Pattern_temp)
        FrequencyArray[j]=FrequencyArray[j]+1
    return FrequencyArray

# #text="AAGCAAAGGTGGG"
# #k=5
# #number=ComputingFrequencies(text,k)
    
def FasterFrequentWords(Text,k):
    FrequentPatterns=[]
    FrequencyArray=ComputingFrequencies(Text,k)
    maxCount=max(FrequencyArray)
    for i in range(0, np.power(4,k)):
        if FrequencyArray[i]==maxCount:
            Pattern=NumberToPattern(i,k)
            FrequentPatterns.append(Pattern)
#    FrequentPatterns= ''.join(FrequentPatterns)
    return FrequentPatterns

# text = "ACTGACTCCCACCCC"
#text = "AAATCGCCCAAACCC"
#tic = time.time()
# pattern=FasterFrequentWords(text, 3)
#toc = time.time()
#shijian = toc-tic
#print(shijian)

def FindingFrequentWorddBySorting(Text,k):
    FrequentPatterns=[]
    Index=[]
    Count=[]
    for i in range(0, (len(Text) - k + 1)):
        Pattern=Text[i : i+k]
        Index.append(PatternToNumber(Pattern))
        Count.append(1)
    SortedIndex=sorted(Index)
    for i in range(0, (len(Text) - k + 1)):
        if SortedIndex[i]==SortedIndex[i-1]:
            Count[i] =Count[i-1] +1
    maxCount=max(Count)
    for i in range(0, (len(Text) - k + 1)):
        if Count[i]==maxCount:
            Pattern=NumberToPattern(SortedIndex[i],k)
            FrequentPatterns.append(Pattern)
    # FrequentPatterns= ''.join(FrequentPatterns)
    return FrequentPatterns

# text="AAGCAAAGGTGGG"
# FrequentPatterns=FindingFrequentWorddBySorting(text,3)

def ClumpFining(Genome,k,t,L):
    FrequentPatterns=[]
    Clump=[]
    for i in range(0, np.power(4,k)):
        Clump.append(0)
    for i in range(0, (len(Genome) - L + 1)):
          Text=Genome[i : i+L] 
          FrequencyArray=ComputingFrequencies(Text,k)
          for index in range(0, np.power(4,k)):
              if FrequencyArray[index]>=t:
                  Clump[index] =1
    for i in range(0, np.power(4,k)):
        if Clump[i]==1:
            Pattern =NumberToPattern(i,k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns

# Genome=data[0:500]
# k=5
# t=3
# L=200
# FrequentPatterns=ClumpFining(Genome,k,t,L)
    
def BetterClumpFining(Genome,k,t,L):
    FrequentPatterns=[]
    Clump=[]
    for i in range(0, np.power(4,k)):
        Clump.append(0)
    Text=Genome[0 : L-1] 
    FrequencyArray=ComputingFrequencies(Text,k)
    for i in range(0, np.power(4,k)):
        if FrequencyArray[i]>=t:
            Clump[i]=1
    for i in range(1, (len(Genome) - L + 1)):
        FirstPattern=Genome[i-1:i-1+k]
        index=PatternToNumber(FirstPattern)
        FrequencyArray[index]=FrequencyArray[index]-1
        LastPattern=Genome[i+L-k-1:i+L-1]
        index=PatternToNumber(LastPattern)
        FrequencyArray[index]=FrequencyArray[index]+1
        if FrequencyArray[index]>=t:
            Clump[index]=1
    for i in range(0, np.power(4,k)):
        if Clump[i]==1:
            Pattern=NumberToPattern(i,k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns

# #Genome=data[0:500]
# #k=5
# #t=3
# #L=200
# #FrequentPatterns=BetterClumpFining(Genome,k,t,L)
   
def Neighbors(Pattern,d):
    if d==0:
        return [Pattern]
    if len(Pattern)==1:
        return ['A','C','G','T']
    Neighborhood=[]
    SuffixNeighbors=Neighbors(Pattern[1:],d)
    for text in SuffixNeighbors:
        if HammingDistance(Pattern[1:], text) < d:
            for n in ['A', 'C', 'G', 'T']:
                Neighborhood.append(n + text)
        else:
            Neighborhood.append(Pattern[0] + text)
    return Neighborhood

# pattern="ACG"
# d=1
# result=Neighbors(pattern, d)

def FrequentWordsWithMismatches(Text,k,d):
    FrequentPatterns=[]
    Close=[]
    FrequencyArray=[]
    for i in range(0, np.power(4,k)):
          Close.append(0)
          FrequencyArray.append(0)
    for i in range(0, (len(Text) - k + 1)):
        Neighborhood=Neighbors(Text[i:i+k],d)
    for Pattern in Neighborhood:
        index=PatternToNumber(Pattern)
        Close[index]=1
    for i in range(0, np.power(4,k)):
        if Close[i] == 1:
            Pattern =NumberToPattern(i,k)
            FrequencyArray[i] =ApproximatePatternCount(Text,Pattern,d)
    maxCount=max(FrequencyArray)
    for i in range(0, np.power(4,k)):
        if FrequencyArray[i] == maxCount:
            Pattern=NumberToPattern(i,k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns

# text="ACACG"
# k=3
# d=1
# pttern=FrequentWordsWithMismatches(text,k,d)

def neighbour(pattern, mismatch, words):
    if mismatch == 0:
        words.add(pattern)
    else:
        bases = ['A', 'T', 'C', 'G']
        for i in range(len(pattern)):
            for j in range(len(bases)):
                new_pattern = pattern[:i] + bases[j] + pattern[i+1:]
                if mismatch <= 1:
                    words.add(new_pattern)
                else:
                    neighbour(new_pattern, mismatch-1, words)

def FrequentWordsWithMismatches_rev(text,k,d):
    allfrequentwords = defaultdict(int)
    for i in range(len(text) - k + 1):
        frequentwords = set()
        neighbour(text[i:i + k], d, frequentwords)
        for words in frequentwords:
            allfrequentwords[words] += 1
    for t in allfrequentwords.keys():
        reverse_k = ReverseComplement(t)
        for i in range(len(text) - k + 1):
            if HammingDistance(text[i:i + k], reverse_k) <= d:
                allfrequentwords[t] += 1
    result = set()
    for t in allfrequentwords.keys():
        if allfrequentwords[t] == max(allfrequentwords.values()):
            result.add(t)
            result.add(ReverseComplement(t))
    for i in result:
        print(i, end=" ")
    return allfrequentwords

# text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
# k, d = 4, 1
# result=FrequentWordsWithMismatches_rev(text,k,d)


def Minimize_Skew(genome):
    skew_value = 0
    min_skew_value = 0
    skew_position = []
    for i, base in enumerate(genome, 1):
        if base.upper() == 'G':
            skew_value = skew_value + 1
        elif base.upper() == 'C':
            skew_value = skew_value - 1
            if skew_value == min_skew_value:
                skew_position.append(i)
            elif skew_value < min_skew_value:
                min_skew_value = skew_value
                skew_position = [i]         
    return skew_position

# position=Minimize_Skew("CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG")

def Find_min_skew(genome):
    skew = [0]
    min_skew = 0
    min_skew_pos = []
    for i, base in enumerate(genome, 1): 
        if base.upper() == 'G':
            skew_val = skew[-1] + 1 
            skew.append(skew_val) 
        elif base.upper() == "C":
            skew_val = skew[-1] - 1
            skew.append(skew_val)
            if skew_val == min_skew:
                min_skew_pos.append(i)
            if skew_val < min_skew:
                min_skew = skew_val 
                min_skew_pos = [i] 
        else:
            skew_val = skew[-1]
            skew.append(skew_val)
    return min_skew_pos, min_skew

# genome = 'CCTATCGGT'
# result=Find_min_skew(genome)

def Draw_GC_picture(genome):
    skew = [0]
    for i, base in enumerate(genome, 1): 
        if base.upper() == 'G':
            skew_val = skew[-1] + 1 
            skew.append(skew_val) 
        elif base.upper() == "C":
            skew_val = skew[-1] - 1
            skew.append(skew_val)
    plt.plot(skew)
    return skew

# genome = data
# result=Draw_GC_picture(genome)

def Compute_prefix_function(text):
    t_len = len(text)
    s = [0] * t_len
    border = 0

    for i in range(1, t_len):
        while (border > 0) and (text[i] != text[border]):
            border = s[border - 1]
        if text[i] == text[border]:
            border += 1
        else:
            border = 0
        s[i] = border
    return s

def Find_pattern(text,pattern):
    """
    Find all the occurrences of the pattern in the text
    and return a list of all positions in the text
    where the pattern starts in the text.
    """
    p_len = len(pattern)
    t_len = len(text)
    result = []
    ref_str = "".join([pattern, "$", text])
    rs_len = len(ref_str)
    prefix_arr = Compute_prefix_function(ref_str)
    for i in range(p_len+1, rs_len):
        if prefix_arr[i] == p_len:
            result.append(i - 2*p_len)
    return result

# text="GATATATGCATATACTT"
# pattern="ATAT"
# result=Find_pattern(text,pattern)

def Find_pattern_Approximate(text,pattern,d):
     return [i\
            for i in range(len(text)-len(pattern)+1)\
            if HammingDistance(pattern,text[i:i+len(pattern)])<=d]

# text="CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC"
# pattern="ATTCTGGA"
# d=3
# result=Find_pattern_Approximate(text,pattern,d)
