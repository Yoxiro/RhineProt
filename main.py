import RhineAMP
LL_37 = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"
aac1 = RhineAMP.PseAAC_Amphiphilic(LL_37,properties=[0,1],Lambda=30,weight=0.05)
# aac = RhineAMP.PseAAC(LL_37,Lambda=30,properties=[0,1])

print(aac1)
