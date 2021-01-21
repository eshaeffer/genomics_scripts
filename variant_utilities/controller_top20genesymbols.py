import argparse
import variant_utilities as var

parser = argparse.ArgumentParser()
parser.add_argument("-input", default="", help="text file containing the paths to the vcf files with each one being on seperate lines")
parser.add_argument("-output", default="", help="Text file containing a path to a vcf file on each line. Also accepts path to the directory containg the vcf files.")
parser.add_argument("-filters", default="", help="Types of filters as an string seperated by commas. Currently accepts the following: passFilter,symbols,supressHeaders")
parser.add_argument("-sort", default="", help="Type of sort as a string. Existing sort types: CaddDescending")
args = parser.parse_args()


if __name__=="__main__":
    loader = var.PickyLoader(args.input)
    err = loader.get_Errors()
    paths = loader.get_fileNames()
    for vcf in paths:
        data = var.PowerfulFormatter(vcf)
        data = var.StrongFilter(data, args.filters).getFormatter()
        var.quickOutput(args.output, vcf, data.getHeaders(), data.getVariants(), args.filters)