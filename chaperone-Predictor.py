import math
import pickle

import numpy
from Bio import SeqIO

with open("./model/std_scale.pkl", 'rb') as file:
    std_scale = pickle.load(file)

with open("./model/model.pkl", 'rb') as file:
    model = pickle.load(file)


class Processor():

    def seqToMat(seq):
        encoder = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                   'Y']
        len = seq.__len__()
        n = int(math.ceil(math.sqrt(len)))
        seqMat = [[0 for x in range(n)] for y in range(n)]
        i = 0
        seqiter = 0
        for i in range(n):
            j = 0
            for j in range(n):
                if seqiter < len:
                    try:
                        aa = int(encoder.index(seq[seqiter]))
                    except ValueError:
                        exit(0)
                    else:
                        seqMat[i][j] = aa
                    seqiter += 1
        return seqMat

    def frequencyVec(seq):
        encoder = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                   'Y']
        fv = [0 for x in range(20)]
        i = 0
        for i in range(20):
            fv[i - 1] = seq.count(encoder[i])
        return fv

    def AAPIV(seq):
        encoder = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                   'Y']
        apv = [0 for x in range(20)]
        i = 1
        sum = 0
        for i in range(20):
            j = 0
            for j in range(len(seq)):
                if seq[j] == encoder[i]:
                    sum = sum + j + 1
            apv[i] = sum
            sum = 0
        return apv[1:] + apv[0:1]

    def print2Dmat(mat):
        n = len(mat)
        i = 0
        strOut = ''
        for i in range(n):
            strOut = strOut + str(mat[i]) + '<br>'
        return strOut

    def PRIM(seq):
        encoder = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                   'Y']
        prim = [[0 for x in range(20)] for y in range(20)]
        i = 0
        for i in range(20):
            aa1 = encoder[i]
            aa1index = -1
            for x in range(len(seq)):
                if seq[x] == aa1:
                    aa1index = x + 1
                    break
            if aa1index != -1:
                j = 0
                for j in range(20):
                    if j != i:
                        aa2 = encoder[j]
                        aa2index = 0
                        for y in range(len(seq)):
                            if seq[y] == aa2:
                                aa2index = aa2index + ((y + 1) - aa1index)
                        prim[i][j] = int(aa2index)
        return prim

    def rawMoments(mat, order):
        n = len(mat)
        rawM = []
        sum = 0
        i = 0
        for i in range(order + 1):
            j = 0
            for j in range(order + 1):
                if i + j <= order:
                    p = 0
                    for p in range(n):
                        q = 0
                        for q in range(n):
                            sum = sum + (((p + 1) ** i) * ((q + 1) ** j) * int(mat[p][q]))
                    rawM.append(sum)
                    sum = 0
        return rawM

    def centralMoments(mat, order, xbar, ybar):
        n = len(mat)
        centM = []
        sum = 0
        i = 0
        for i in range(order + 1):
            j = 0
            for j in range(order + 1):
                if i + j <= order:
                    p = 0
                    for p in range(n):
                        q = 0
                        for q in range(n):
                            sum = sum + ((((p + 1) - xbar) ** i) * (((q + 1) - ybar) ** j) * mat[p][q])
                    centM.append(sum)
                    sum = 0
        return centM

    def hahnMoments(mat, order):
        N = len(mat)
        hahnM = []
        i = 0
        for i in range(order + 1):
            j = 0
            for j in range(order + 1):
                if i + j <= order:
                    answer = Processor.hahnMoment(i, j, N, mat)
                    hahnM.append(answer)
        return hahnM

    def hahnMoment(m, n, N, mat):
        value = 0.0
        x = 0
        for x in range(N):
            y = 0
            for y in range(N):
                value = value + (
                        mat[x][y] * (Processor.hahnProcessor(x, m, N)) * (Processor.hahnProcessor(x, n, N)))
        return value

    def hahnProcessor(x, n, N):
        return Processor.hahnPol(x, n, N) * math.sqrt(Processor.roho(x, n, N))

    def hahnPol(x, n, N):
        answer = 0.0
        ans1 = Processor.pochHammer(N - 1.0, n) * Processor.pochHammer(N - 1.0, n)
        ans2 = 0.0
        k = 0
        for k in range(n + 1):
            ans2 = ans2 + math.pow(-1.0, k) * ((Processor.pochHammer(-n, k) * Processor.pochHammer(-x, k) *
                                                Processor.pochHammer(2 * N - n - 1.0, k)))
        answer = ans1 + ans2
        return answer

    def roho(x, n, N):
        return Processor.gamma(n + 1.0) * Processor.gamma(n + 1.0) * Processor.pochHammer((n + 1.0),
                                                                                          N)

    def gamma(x):
        return math.exp(Processor.logGamma(x))

    def logGamma(x):
        temp = (x - 0.5) * math.log(x + 4.5) - (x + 4.5)
        ser = 101.19539853003
        return temp + math.log(ser * math.sqrt(2 * math.pi))

    def pochHammer(a, k):
        answer = 1.0
        i = 0
        for i in range(k):
            answer = answer * (a + i)
        return answer

    def calcFV(seq):
        fv = [0 for x in range(150)]
        fvIter = 0
        myMat = Processor.seqToMat(seq)
        myRawMoments = Processor.rawMoments(myMat, 3)
        for ele in myRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        xbar = myRawMoments[4]
        ybar = myRawMoments[1]
        myCentralMoments = Processor.centralMoments(myMat, 3, xbar, ybar)
        for ele in myCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        myHahnMoments = Processor.hahnMoments(myMat, 3)
        for ele in myHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        myFrequencyVec = Processor.frequencyVec(seq)
        for ele in myFrequencyVec:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        myPRIM = Processor.PRIM(seq)
        myPRIMRawMoments = Processor.rawMoments(myPRIM, 3)
        xbar2 = myPRIMRawMoments[4]
        ybar2 = myPRIMRawMoments[1]
        myPRIMCentralMoments = Processor.centralMoments(myPRIM, 3, xbar2, ybar2)
        for ele in myPRIMRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        for ele in myPRIMCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        myPRIMHahnMoments = Processor.hahnMoments(myPRIM, 3)
        for ele in myPRIMHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        myAAPIV = Processor.AAPIV(seq)
        for ele in myAAPIV:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        myRPRIM = Processor.PRIM(seq[::-1])
        myRPRIMRawMoments = Processor.rawMoments(myRPRIM, 3)
        xbar3 = myRPRIMRawMoments[4]
        ybar3 = myRPRIMRawMoments[1]
        myRPRIMCentralMoments = Processor.centralMoments(myRPRIM, 3, xbar3, ybar3)
        for ele in myRPRIMRawMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        for ele in myRPRIMCentralMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        myRPRIMHahnMoments = Processor.hahnMoments(myRPRIM, 3)
        for ele in myRPRIMHahnMoments:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        myRAAPIV = Processor.AAPIV(seq[::-1])
        for ele in myRAAPIV:
            fv[fvIter] = ele
            fvIter = fvIter + 1
        fv = numpy.asarray(fv).reshape(1, -1)
        return fv

    def matrix(string, length):
        return (string[0 + i:length + i] for i in range(0, len(string), length))

    def driverMethod(seq):
        if seq.isalpha():
            out = str(Processor.performPrediction(numpy.asarray(Processor.calcFV(seq))))
            seq = "\n".join(list(Processor.matrix(seq, 100))) + '\n'
            return [seq, out]
        else:
            return ['Invalid Sequence']

    def scaling(FV):
        FV = std_scale.transform(FV)
        return FV

    def performPrediction(FV):
        FV = Processor.scaling(FV)
        output = model.predict(FV)
        if output == 1:
            return 'Chaperone'
        else:
            return 'Non-Chaperone'


i = 0
for seq_record in SeqIO.parse("./data/chaperone2980.fasta", "fasta"):
    seq = str(seq_record.seq)
    print(str(i) + ': ' + str(Processor.driverMethod(seq)))
    i = i + 1
