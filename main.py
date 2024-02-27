import RhineAMP
LL_37 = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"
aac = RhineAMP.PseAAC(LL_37,Lambda=30,weight=0.05)
aac1 = RhineAMP.PseAAC_Amphiphilic(LL_37,Lambda=30,weight=0.05)

print(aac1)
