#!/usr/bin/env python
import sys
import argparse
import os

def main():
    parser = argparse.ArgumentParser(prog='globalalign',description='find the global alignment of two sequences')
    parser.add_argument('sequence1', type = str, help = 'the first sequence') #file location
    parser.add_argument('sequence2', type = str, help = 'the second sequence')
    parser.add_argument('--match', type = int, default = 1)
    parser.add_argument('--mismatch', type = int, default = -1)
    parser.add_argument('--gap', type = int, default = -1, help = 'Gap penalty')
    parser.add_argument('--similarity_matrix', type = str)

    args = parser.parse_args()
    file1 = args.sequence1
    file2 = args.sequence2
    match = args.match
    mismatch = args.mismatch
    indel = args.gap
    similarity_matrix = args.similarity_matrix

    # need to add in check if file exists
    matrix = {}
    if(similarity_matrix is not None):
        dirname = os.path.dirname(__file__)
        filenameS = os.path.join(dirname, 'matrices/',similarity_matrix)
        with open(filenameS, 'r') as reader:
            lines = reader.readlines()
            for l in lines:
                if(l[0] == "#"):
                    continue
                elif(l[0] == " "):
                    ls = l.split()
                    # print(ls)
                    continue
                else:
                    vals = l.split()
                    if(len(vals) == 0):
                        continue
                    letter = vals[0]
                    for index, item in enumerate(ls):
                        matrix[(letter,item)] = vals[index+1]
 
    def similarityMatrix(a,b): #basic similarity matrix
        # print(a,b)
        if (a=='U'):
            a = 'T'
        if (b=='U'):
            b = 'T'
        if(similarity_matrix is None):
            if(a==b):
                return match
            else:
                return mismatch
        else:
            if(a == " " or b == " "):
                return mismatch
            return int(matrix.get((a,b)))


    with open(file1, 'r') as reader:
        sequence1 = reader.readline()
        sequenceLength1 = len(sequence1)
        if(sequence1[len(sequence1)-1] == '\n'):
            sequenceLength1 = sequenceLength1 - 1
            sequence1 = sequence1[:sequenceLength1]
        print("Sequence 1")
        print(sequence1)

    with open(file2, 'r') as reader:
        sequence2 = reader.readline()
        sequenceLength2 = len(sequence2)
        if(sequence2[len(sequence2)-1] == '\n'):
            sequenceLength2 = sequenceLength2 - 1
            sequence2 = sequence2[:sequenceLength2]
        print("Sequence 2")
        print(sequence2)

    rows, cols = (sequenceLength2 + 1, sequenceLength1 + 1)
    Fmatrix = [[0 for i in range(cols)] for j in range(rows)]
    for i in range(0, sequenceLength2 + 1):
        Fmatrix[i][0] = indel * i
    for j in range(0, sequenceLength1 + 1):
        Fmatrix[0][j] = indel * j
    for i in range(1, sequenceLength2 + 1):
        for j in range(1, sequenceLength1 + 1):
                M = Fmatrix[i-1][j-1] + similarityMatrix(sequence2[i-1], sequence1[j-1])
                D = Fmatrix[i-1][j] + indel
                I = Fmatrix[i][j-1] + indel
                Fmatrix[i][j] = max(M,D,I)

    AllignmentA = ""
    AllignmentB = ""
    global allCount
    allCount = 0

    def findAllignment(AllA,AllB,i_x,j_x,score,output_file):
        
        if (i_x == 0 and j_x == 0):
            global allCount
            allCount = allCount + 1
            print("Best Allignment", allCount)
            print(AllB)
            print(AllA)
            output_file.writelines(["Best Allignment ", str(allCount), "\n", AllB, "\n", AllA, "\n"])
            return
        nextScores = [Fmatrix[i_x-1][j_x] + indel, Fmatrix[i_x][j_x-1] + indel, Fmatrix[i_x-1][j_x-1]+ similarityMatrix(sequence2[i-1],sequence1[j-1])]
        mNext = max(nextScores)
        newScore = mNext + score

        #if there is a match or mismatch
        if (i_x> 0 and j_x > 0 and Fmatrix[i_x][j_x] == Fmatrix[i_x-1][j_x-1] + similarityMatrix(sequence2[i_x-1],sequence1[j_x-1])):
            AllAn = sequence2[i_x - 1] + AllA
            AllBn = sequence1[j_x - 1] + AllB
            findAllignment(AllAn,AllBn,i_x-1,j_x-1,newScore,output_file)
        #if there is a deletion then match allignmentA with a gap
        if (i_x > 0 and Fmatrix[i_x][j_x] == Fmatrix[i_x-1][j_x] + indel):
            AllAn = sequence2[i_x - 1] + AllA
            AllBn = "-" + AllB
            findAllignment(AllAn,AllBn,i_x-1,j_x,newScore,output_file)
        #if there is a insertion then match allignmentB with a gap
        if (j_x > 0 and Fmatrix[i_x][j_x] == Fmatrix[i_x][j_x-1] + indel):
            AllAn = "-" + AllA
            AllBn = sequence1[j_x - 1] + AllB
            findAllignment(AllAn,AllBn,i_x,j_x-1,newScore,output_file)

    with open ("Output.txt", "w") as output_file:
        output_file.writelines(["Sequence 1","\n",sequence1,"\n","Sequence 2","\n",sequence2,"\n"])
        findAllignment("","",sequenceLength2,sequenceLength1,0,output_file)
    # print("Global Allignment")
    # print(AllignmentB)
    # print(AllignmentA)
    # with open ("Output.txt", "w") as output_file:
    #     output_file.writelines([AllignmentB,"\n", AllignmentA])
if __name__ == '__main__':
    main()
