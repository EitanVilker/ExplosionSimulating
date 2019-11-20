import maya.cmds as cmds

filePath = "C:/Users/naman/Downloads/ExplosionSimulating/NURBS/basic.txt"
fileHandler = open(filePath, "w")    

n = 200
for x in range(n):
    c = float(x) / n
    fileHandler.write(",".join(map(str,cmds.pointOnCurve('curve1', pr=c, top=True))))
    fileHandler.write("\n")
fileHandler.close()