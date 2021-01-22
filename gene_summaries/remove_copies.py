import argparse

#
# Argparse to choose arguments
#
parser = argparse.ArgumentParser()
parser.add_argument('-geneList', default='genes.txt', help='file consisting of a single column of gene names with a column header')
parser.add_argument('-geneSummaries', default='gene_summaries.txt', help='file consisting of 4 columns seperate by tabs - first being a number, second the gene name, and 3 and 4 are summaries')
parser.add_argument('-output', default='summaries_noRepeats.txt', help='what you want your output to be called - output consists of gene summaries, but it matches your gene list')
args = parser.parse_args()


#
# Function to get the lines of a file and return as list
#
def getFileLines(filePath):
   f = open(filePath, 'r')
   line = f.readline()
   lineList = []
   while(line != ""):
      lineList.append(line)
      line = f.readline()
   f.close()
   return lineList

#
# Function to remove consecutive repeats from second argument
# until it matches with first argument. Returns non-repeated list
#
def deleteConsecutiveRepeats(lines, fixedList):
   for i in range(len(lines)):
      # Ignore first line since it contains the headers
      if(i == 0):
         pass
      else:
         # Here we are getting the thing we want to match
         line = lines[i].replace("\n","")
         try:
            lineToCheck = fixedList[i].split('\t')[1].replace('\"',"")
         except:
            pass
         while(line != lineToCheck):
            try:
               del fixedList[i]
               # Again, will likely need to be changed to accomidate different files
               lineToCheck = fixedList[i].split('\t')[1].replace('\"',"")
            except:
               break
   return fixedList

#
# Outputs a list to a filename of choice
#
def outputLines(lines, outputName):
    f = open(outputName, 'w+')
    for line in lines:
        f.write(line)
    f.close()


if __name__ == "__main__":
   genes = getFileLines(args.geneList)
   summaries = getFileLines(args.geneSummaries)
   fixedSummaries = deleteConsecutiveRepeats(genes, summaries)
   outputLines(fixedSummaries, args.output)
