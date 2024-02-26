import RhineAMP
LL_37 = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"

from RhineAMP.Descriptors.PseAAC import normalized_aac

nor = normalized_aac(LL_37)
std = nor.std(ddof=0)
mean = nor.mean()
print(nor)
print(std)
print(mean)