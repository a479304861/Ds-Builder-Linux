import os
class processDrimm:
    def __init__(self, inputSequencePath, output,cycleLength,dustLength,speciesAndChrLen):
      self.inputSequencePath = inputSequencePath
      self.output = output
      self.cycleLength = cycleLength
      self.dustLength = dustLength
      self.speciesAndChrLen = speciesAndChrLen

    def execute(self):
        AbsSequencePath=os.path.abspath(self.inputSequencePath)
        AbsSequencePath= AbsSequencePath.replace("\\", "/")
        AbsOutPath=os.path.abspath(self.output)
        AbsOutPath=AbsOutPath.replace("\\", "/")
        cmd = "mono ./utils/Program.exe "+AbsSequencePath+" "+AbsOutPath+" "+self.cycleLength+" "+self.dustLength
        # main = "mono ./win-x64/drimm.dll "+AbsSequencePath+" "+AbsOutPath+" "+self.cycleLength+" "+self.dustLength    linux
        print(cmd)
        drimmPrint = os.system(cmd)
        print (drimmPrint)
        self.splitBySpecies()

    def splitBySpecies(self):
        readfile = open(self.output+"/blocks.txt",'r')
        blocks= readfile.read().split('\n')
        readfile.close()
        speciesBlock={}
        index=0
        for speciesKey in self.speciesAndChrLen:
            speciesBlock[speciesKey]=blocks[index:index+self.speciesAndChrLen[speciesKey]]
            index+=self.speciesAndChrLen[speciesKey]
        for speciesBlockItem in speciesBlock:
            output = open(self.output+"/"+speciesBlockItem+".block",'w')
            for s in speciesBlock[speciesBlockItem]:
                output.write(s)
                output.write("\n")
            output.close()
        


