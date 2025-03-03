import sys     
import numpy as np
import matplotlib.pyplot as plt

def readData(fname):                              #Defined a function to read data from given data set
    print(fname)
    with open(fname, 'r') as f:
        lines = f.readlines()
    V = float(lines[0].split(':')[-1].strip())

    Sigma = []
    for line in lines[1:]:
        
        m_list = []
        for sigma in line.split():
            a = float(sigma.strip())
            m_list.append(a)

        Sigma.append(m_list)
            

    Sigma = np.array(Sigma)
    return Sigma, V

def AvgStress(n,SigArr,VolumeArr):
    result = 0
    for i in range(n):
        result = result + SigArr[i] * VolumeArr[i]
    
    return result  / sum(VolumeArr[:n])

def EigVal(Sigma):                                    #Eigen values function 
    return np.linalg.eigvals(Sigma)

def plotMohrsCircle(SigmaEigval, foutNamePNG):        #Ploting Mohr's Circle
    print(foutNamePNG)

    # convert MPa to GPa
    SigmaEigval = list(reversed(sorted([ v / 1e3 for v in SigmaEigval])))

    sig1 = (SigmaEigval[1] + SigmaEigval[2]) / 2
    sig2 = (SigmaEigval[0] + SigmaEigval[2]) / 2
    sig3 = (SigmaEigval[0] + SigmaEigval[1]) / 2

    R1 = (SigmaEigval[1] - SigmaEigval[2]) / 2
    R2 = (SigmaEigval[0] - SigmaEigval[2]) / 2
    R3 = (SigmaEigval[0] - SigmaEigval[1]) / 2

    print(sig1, R1)
    print(sig2, R2)
    print(sig3, R3)
    

    C2 = plt.Circle((sig2, 0), R2, facecolor='#135198', edgecolor='black', antialiased=True, linewidth=1.0)
    C1 = plt.Circle((sig1, 0), R1, facecolor='white', edgecolor='black', antialiased=True, linewidth=1.0)
    C3 = plt.Circle((sig3, 0), R3, facecolor='white', edgecolor='black', antialiased=True, linewidth=1.0)

    fig, ax = plt.subplots(figsize=(6.4, 6.4))

    ax.add_patch(C2)
    ax.add_patch(C1)
    ax.add_patch(C3)

    ax.set_ylabel(r'$\tau$ [GPa]')
    ax.set_xlabel(r'$\sigma$ [GPa]')
    

    ax.autoscale_view()

    fig.savefig(foutNamePNG, pad_inches=0.5)

def main(FilePrefixStr, MatPtStr, FileExtension, NumMaterialPoints):

    SigmaArr = []
    VolumeArr = []
    for n in range(1, NumMaterialPoints + 1):
        Sigma, V = readData(f"{FilePrefixStr}_{MatPtStr}{n}.{FileExtension}")   #Using f strings
        SigmaArr.append(Sigma)
        VolumeArr.append(V)

    AvgSigma = AvgStress(NumMaterialPoints, SigmaArr, [1.0] * NumMaterialPoints)
    VolAvgSigma = AvgStress(NumMaterialPoints, SigmaArr, VolumeArr)
    
    for n in range(1, NumMaterialPoints + 1):
        SigmaEigval = EigVal(SigmaArr[n - 1])

        plotMohrsCircle(SigmaEigval, f"MohrCirc_{FilePrefixStr}_{MatPtStr}{n}.PNG")

    AvgSigmaEigval = EigVal(AvgSigma)
    plotMohrsCircle(AvgSigmaEigval, f"MohrCirc_{FilePrefixStr}_Avg.PNG")

    VolAvgSigmaEigval = EigVal(VolAvgSigma)
    plotMohrsCircle(VolAvgSigmaEigval, f"MohrCirc_{FilePrefixStr}_VolAvg.PNG")


if __name__ == '__main__':
    
    if len(sys.argv) != 5:
        print("Usage: python MohrCirc.py <FilePrefixStr> <MatPtStr> <FileExtension> <NumMaterialPoints>")
        exit(1)
    
    FilePrefixStr = sys.argv[1]
    MatPtStr = sys.argv[2]
    FileExtension = sys.argv[3]
    NumMaterialPoints = int(sys.argv[4])
    main(FilePrefixStr, MatPtStr, FileExtension, NumMaterialPoints)
