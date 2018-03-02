# define biological parameters used for simulation
# left colon

# Males
slope = 10.42e-5 #7.34e-5
X = 1.e8
qI = 1.52e-5 #3.66e-5

# net cell proliferation
gI = 0.144 # 0.122 #; gM = 2.4
pI = -gI

# sym cell division
alphaI = 9; pinfI = gI/alphaI
# alphaM = 100; pinfM = gM/alphaM
gamI = 1-pinfI 

mu0 = mu1 = sqrt(slope/X/pinfI)
# mu2eff = qI*pinfI*pinfM; mu2 = mu2eff/pinfM
# pM = -gM; qM = rho/pinfM
# tlag = -log(alphaM*rho/gM/gM)/gM

tlag=0



