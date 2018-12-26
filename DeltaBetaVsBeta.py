import numpy as np 
import matplotlib.pyplot as plt 
plt.switch_backend('Agg')
import glob


files = glob.glob('DiagFiles/*Beta*SameK100Sims.txt')

TrueBetas = []
DeltaBetas = []
SigmaBetas = []

files = sorted(files)

for f in files:
	if 'MCBeta' in f:
		continue
	BetaFile = np.genfromtxt(f, names=True, dtype = None, skip_header = 10)
	print BetaFile.shape
	try:
		TrueBetas.append(BetaFile['DataBeta'][-1])
		DeltaBetas.append(BetaFile['delta_Beta'][-1])
		SigmaBetas.append(BetaFile['sigma_Beta'][-1])
	except:
		TrueBetas.append(float(BetaFile['DataBeta']))
		DeltaBetas.append(float(BetaFile['delta_Beta']))
		SigmaBetas.append(float(BetaFile['sigma_Beta']))


print TrueBetas
print DeltaBetas

plt.xlabel("True Data Beta")
plt.ylabel("Delta Beta")
plt.scatter(TrueBetas, np.array(DeltaBetas) + 2.11 - np.array(TrueBetas), marker = 'o')
plt.errorbar(TrueBetas, np.array(DeltaBetas) + 2.11 - np.array(TrueBetas), yerr = SigmaBetas, fmt = 'none', marker = 'o')
#plt.plot(TrueBetas, np.array(TrueBetas) - 2.11)
plt.plot(TrueBetas, np.zeros(len(TrueBetas)) )
plt.savefig('DeltaBetaVsBeta_NoCheat.png')

plt.clf()

files = glob.glob('DiagFiles/*Beta*SameK100Sims.txt')

TrueBetas = []
DeltaBetas = []
SigmaBetas = []

files = sorted(files)

for f in files:
	if not ('MCBeta15' in f):
		continue
	print f
	BetaFile = np.genfromtxt(f, names=True, dtype = None, skip_header = 10)
	print BetaFile.shape
	try:
		TrueBetas.append(BetaFile['DataBeta'][-1])
		DeltaBetas.append(BetaFile['delta_Beta'][-1])
		SigmaBetas.append(BetaFile['sigma_Beta'][-1])
	except:
		TrueBetas.append(float(BetaFile['DataBeta']))
		DeltaBetas.append(float(BetaFile['delta_Beta']))
		SigmaBetas.append(float(BetaFile['sigma_Beta']))


print TrueBetas
print DeltaBetas

plt.xlabel("True Data Beta")
plt.ylabel("Delta Beta")
plt.scatter(TrueBetas, np.array(DeltaBetas) + 1.51 - np.array(TrueBetas), marker = 'o')
plt.errorbar(TrueBetas, np.array(DeltaBetas) + 1.51 - np.array(TrueBetas), yerr = SigmaBetas, fmt = 'none', marker = 'o')
#plt.plot(TrueBetas, np.array(TrueBetas) - 1.51)
plt.plot(TrueBetas, np.zeros(len(TrueBetas)) )
plt.savefig('DeltaBetaVsBeta_MCBeta15_NoCheat.png')