import os, re
import numpy as np
import subprocess
import traceback as trace
import pandas as pd

# class PickyLoader, used to get filenames and get sample names from filenames, requires path to a directory/file to
# instantiate object has private member fileNames which are the fileNames inputted as arguments, sampleNames which
# are all the sampleNames, and errorLog which holds the traces for each method has getters/setters for above instance
# variables as well as a methods for retrieving the filenames: loadTheFile, loadTheDir. It also has method for
# getting filenames: getSampleNames
class PickyLoader:
   def __init__(self, inpt):
      self.__fileNames = ""
      self.__sampleNames = []
      self.__errorLog = []
      if os.path.isdir(inpt):
         self.loadTheDir(inpt)
      else:
         self.loadTheFile(inpt)
      self.getSampleNames()

   def get_fileNames(self):
      return self.__fileNames
   def set_fileNames(self, inpt):
      self.__fileNames = inpt
   def get_sampleNames(self):
      return self.__sampleNames
   def set_sampleNames(self, inpt):
      self.__sampleNames = inpt
   def get_Errors(self):
      return self.__errorLog
   def add_error(self,error):
      self.__errorLog.append(error)

   # loadTheDir: a way to get all root paths to each file in a directory
   # Expects the path to a valid directory name as its input
   # Has no outputs - it simply sets the fileNames variable
   def loadTheDir(self, inpt):
      tmp = []
      vcf = re.compile('.vcf$')
      try:
         for root, dirs, files in os.walk(inpt):
            for file in files:
               m = vcf.search(file)
               if m:
                  tmp.append(os.path.join(root,file))
               else:
                  print("Ommitting the following non-vcf file found in inputted directory: " + file)
         tmp.sort()
         self.__fileNames = np.array(tmp)
         if np.size(self.__fileNames) == 0:
            print("warning - no vcfs detected in the inputted directory: " + inpt)
      except:
         print("Something went wrong getting the file-paths from the directory")
         self.__errorLog.append(trace.format_exc())

   # loadTheFile: a simple way to create an array of filenames using a text file with files on each line
   # Expects as an input the path to a text file containing the complete filenames
   # Has no outputs - appends directly to the self.__fileNames instance variable
   def loadTheFile(self, inpt):
      tmp = []
      try:
         tmp = np.loadtxt(inpt, dtype="str")

         for files in np.nditer(tmp):
            # need to be string
            files = np.array2string(files).strip().replace('\'','')
            if not os.path.isfile(files):
               print("Ommiting the following non-existent file -- " + files)
               tmp = np.delete(tmp, np.where(tmp == files))
         tmp = np.sort(tmp)
         self.__fileNames = tmp
      except:
         self.__fileNames = np.array([""])
         print("something went wrong getting the file-paths from the text file.")
         self.__errorLog.append(trace.format_exc())

   # getSampleNames: Creates the sample-names from the fileNames instance variable
   # It does this by using basename() to get the filename, then splitting the filename,
   # taking the first element. As such it is pretty picky, hence the name of the class.
   def getSampleNames(self):
      tmp = []
      for path in np.nditer(self.__fileNames):
         try:
            path = np.array2string(path).replace('\'','')
            tmp1 = os.path.basename(path).split(".")[0]
            tmp.append(tmp1)
         except:
            print("something went wrong getting the sample name of: " + path)
            print("file names must be in the following format: SAMPLE-NAME.OPTIONS.vcf")
            self.__errorLog.append(trace.format_exc())
      self.__sampleNames = np.array(tmp)


# PowerfuleFormatter class: stores headers and variants, variants in a pandas dataframe requires types of filters and
# the path to the file Has private instance variants path which is the path to the vcf file to filter, variants which
# is a pandas dataframe containing filtered/raw variants headers which is the headers to the vcf file to filter,
# and types which is all the filters to run on this particular vcf file.
class PowerfulFormatter:
   def __init__(self, path, numpy = False):
      self.__variants = ""
      self.__headers = ""
      self.dataGrabber(path, numpy)
      self.__path = path
      self.__errorLog = []

   def get_errors(self):
      return self.__errorLog
   def getVariants(self):
      return self.__variants
   def setVariants(self, variants):
      self.__variants = variants
   def getHeaders(self):
      return self.__headers
   def setHeaders(self, headers):
      self.__headers = headers
   def add_error(self, error):
      self.__errorLog.append(error)
   def get_path(self):
      return self.__path

   # Data Grabber takes path to a vcf file as an input.
   # Uses shell commands to get info about vcf (headers, variants)
   # invokes the frameData put data in pandas dataframe with 10 columns and N rows where N is the number of variants.
   def dataGrabber(self, path, numpy):
      headers = "cat " + path + " | grep \"##\" "
      header = "cat " + path + " | grep \"#CHROM\" "
      variants = "cat " + path + " | grep \"^chr\" "
      try:
         self.__headers = np.array(subprocess.check_output(headers, shell=True).decode("utf-8").split("\n"))
         self.__headers.resize(np.size(self.__headers) - 1)
         self.__variants = subprocess.check_output(variants, shell=True).decode("utf-8").split("\n")
         data = subprocess.check_output(variants, shell=True).decode("utf-8")
         header = subprocess.check_output(header, shell=True).decode("utf-8")
         if numpy:
            data = data.split("\n")
            del data[-1]
            header = header.rstrip()
            header = header.split("\t")
            tmp = []
            tmp.append(header)
            for i in data:
               tmp.append(i.split("\t"))
            self.__variants = np.array(tmp)
         else:
            self.frameData(header)

      except:
         print("Problem using system calls to collect headers/variants in the \n"
               "following location: " + path)
         self.__errorLog.append(trace.format_exc())

   # Data framer creates a pandas data frame based on what dataGrabber found in the vcf file
   # It takes in the header to the variants and then creates a pandas dataframe based on the
   # variants dataGrabber had found.
   def frameData(self, header):
      count = 0
      del self.__variants[-1]
      for i in self.__variants:
         self.__variants[count] = i.split("\t")
         count += 1
      self.__variants = np.array(self.__variants)
      header = header.split("\t")
      header[-1] = header[-1].replace("\n","")
      try:
         self.__variants = pd.DataFrame(self.__variants, columns=header)
      except:
         print("An issue occured when framing data for the following sample: " + header[-1])
         self.__errorLog.append(trace.format_exc())


# Class Strong filter is used to filter variants based on information located in the vcf file
# Takes in a formatter class and a list of filters to apply
# Has one filter currently, passFilter
class StrongFilter:
   def __init__(self, formatter, types=""):
      filterDictionary = {'passFilter': self.passFilter, 'symbols': self.symbols, 'supressHeaders': self.supressHeaders}
      self.__formatter = formatter
      self.__errorLog = []
      types = types.split(",")
      if types[0] == "":
         print("No filters provided for instantiated filters Object, exiting without filtering.")
      else:
         for i in types:
            filterDictionary[i]()

   def getFormatter(self):
      return self.__formatter
   def getErrorLog(self):
      return self.__errorLog

   #passFilter creates sets the variants to only those that are passed
   def passFilter(self):
      vars = self.__formatter.getVariants()
      try:
         self.__formatter.setVariants(vars[vars['FILTER'] == 'PASS'])
         self.__formatter.setHeaders(np.append(self.__formatter.getHeaders(),"##Pass Filtered by variant_utilites"))
      except:
         print("error occured while pass-Filtering: " + list(vars.columns.values)[-1])
         self.__errorLog.append(trace.format_exc())

   def symbols(self):
      vars = self.__formatter.getVariants()
      try:
         vars = vars['INFO']
         vars = pd.DataFrame(vars)
         vars = vars.applymap(lambda x: x.split("|"))
         vars = vars.applymap(lambda x: x[3])
         self.__formatter.setVariants(vars)
         self.__formatter.setHeaders(np.append(self.__formatter.getHeaders,"##Gene symbols only by variant_utilities"))
      except:
         print("error occured while filtering by gene symbols: ")
         self.__formatter.add_error(trace.format_exc())

   def supressHeaders(self):
      self.__formatter.setHeaders(np.array([]))



# Class StrangeSorter exists to provide support genomics analysis by sorting vcf files according to certain metrics.
# Currently existing sort types: by descending CADD score. Any other types must be added to dict to work!
# outputs nothing, but mutates the inputted formatter to be used in further analysis/ written to a file.
class StrangeSorter:
   def __init__(self,formatter, type = ""):
      sortDictionary = {'CaddDescending': self.byCaddDescending}
      self.__formatter = formatter
      self.__errorLog = []

      if type == "":
         print("No types provided, StrangeSorter cannot provide sorted files w/out type of sort.")
      else:
         try:
            sortDictionary[type]()
         except:
            print("The following sort order does not exist: " + type)
            self.__errorLog.append(trace.format_exc())

   def getFormatter(self):
      return self.__formatter
   def setFormatter(self,formatter):
      self.__formatter = formatter

   def byCaddDescending(self):
      newFrame = pd.DataFrame({"CADD": self.__formatter.getVariants()['INFO']})
      newFrame = newFrame.applymap(lambda x: x.split('|'))
      newFrame = newFrame.applymap(lambda x: x[len(x) - 5])
      self.__formatter.setVariants(self.__formatter.getVariants().join(newFrame))
      self.__formatter.setVariants(self.__formatter.getVariants().sort_values(by=['CADD'],ascending=False).drop(['CADD'], axis=1))
      self.__formatter.setHeaders(np.append(self.__formatter.getHeaders(),"##Sorted by descending CADD score using variant utilities"))

# This is an example of an idiomatic python counting algorithm. It uses the information encoded in a formatter class
# to allow you to output stats about the variants directly to a xls format. To instantiate one, you need to provide
# a formatter with it's variants and headers members filled. It also needs the 'type' you'd like to count along with
# their cutoffs (currently accepts cadd, maf, noMaf, noCadd). recommended cutoffs for each: 20, 0.01. You need to provide a filename to
# output to and sample's name as well.
class ZenCounter:
   def __init__(self, formatter, typesToCount, cutoffs, filename, sampleName):
      countDict = {'cadd': self.cadd, 'maf': self.maf}
      self.__formatter = formatter
      self.__types = typesToCount
      self.__cutoff = cutoffs
      self.__filename = filename
      self.__sampleName = sampleName


#####################################
############ Helpers ################
#####################################

# Helper function limit returns the top n number of entries in the variant call format files
# Takes in the number of entries to return as an integer, and the formatter in which the variants are contained
# returns the formatted formatter class with only n entries in the variants member.
def limit(data, n):
   try:
      data.setVariants(data.getVariants().head(n))
      data.setHeaders(np.append(data.getHeaders(), "##Output limited via variant utilities"))
   except:
      print("Problem limiting output")
      data.add_error(trace.format_exc())
   return data

# quickOutput writes a file using headers and variant from previous functions
# takes in outputPath, inputPath with fileName, headers, variants, and filters used
# outputs nothing, but creates a file in outPath, with the same name as it, but augmented with filters used.
def quickOutput(outPath, inPath, headers, variants, filters="", sort=""):

   name = os.path.basename(inPath).replace("vcf","") + ".".join(filters.split(",")) + sort + ".vcf"
   outPath += "/" + name
   remove = "rm variants1.txt headers.txt variants.txt"
   cut = "cut -f2- variants1.txt > variants.txt"
   concatFiles = "cat headers.txt variants.txt >> " + outPath
   try:
      np.savetxt("headers.txt",headers, fmt="%s")
   except:
      print("error outputing headers from: " + name)
      print(trace.format_exc())
   try:
      variants.to_csv('variants1.txt', sep="\t")
   except:
      print("error outputting vcfs from: " + name)
      print(trace.format_exc())

   subprocess.check_output(cut, shell=True)
   subprocess.check_output(concatFiles, shell=True)
   subprocess.check_output(remove, shell=True)