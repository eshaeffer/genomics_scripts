import pandas as pd
import numpy as np
import variant_utilities as var

if __name__=="__main__":
    summary = pd.read_table("/Users/eliasshaeffer/Desktop/somatic_gatk_low/summaries_noRepeats.txt", sep='\t')
    loader = var.PickyLoader('/Users/eliasshaeffer/Desktop/somatic_gatk_low')
    sampleNames = loader.get_sampleNames()
    files = loader.get_fileNames()
    data = var.PowerfulFormatter(files[0])
    data = var.limit(data, 20)

    newFrame = pd.DataFrame({"CADD": data.getVariants()['INFO']})
    newFrame = newFrame.applymap(lambda x: x.split('|'))
    newFrame = newFrame.applymap(lambda x: x[len(x) - 6])
    data.setVariants(data.getVariants().join(newFrame))

    newFrame2 = pd.DataFrame({"MAF": data.getVariants()['INFO']})
    newFrame2 = newFrame2.applymap(lambda x: x.split(';'))
    newFrame2 = newFrame2.applymap(lambda x: x[-2].replace("bcfGnomAD_AF=",""))
    data.setVariants(data.getVariants().join(newFrame2))

    sampleArray = np.repeat(sampleNames[0], 20)
    sampleFrame = pd.DataFrame(sampleArray, columns=['SampleNames'])
    data.setVariants(data.getVariants().join(sampleFrame))

    vars = data.getVariants()#.drop(data.getVariants().columns[len(data.getVariants().columns) - 3], axis=1)


    sampleNames = np.delete(sampleNames, 0)
    files = np.delete(files, 0)
    for i in range(0, np.size(files)):
        data = var.PowerfulFormatter(files[i])
        newFrame = pd.DataFrame({"CADD": data.getVariants()['INFO']})
        newFrame = newFrame.applymap(lambda x: x.split('|'))
        newFrame = newFrame.applymap(lambda x: x[len(x) - 6])
        data.setVariants(data.getVariants().join(newFrame))

        newFrame2 = pd.DataFrame({"MAF": data.getVariants()['INFO']})
        newFrame2 = newFrame2.applymap(lambda x: x.split(';'))
        newFrame2 = newFrame2.applymap(lambda x: x[-2].replace("bcfGnomAD_AF=", ""))
        data.setVariants(data.getVariants().join(newFrame2))


        sampleArray = np.repeat(sampleNames[i], 20)
        sampleFrame = pd.DataFrame(sampleArray, columns=['SampleNames'])
        data.setVariants(data.getVariants().join(sampleFrame))

        vars = pd.DataFrame.append(vars, data.getVariants())#.drop(data.getVariants().columns[len(data.getVariants().columns) - 1], axis=1))
    vars = pd.DataFrame.reset_index(vars)
    combined = pd.concat([summary, vars], axis=1)
    combined.to_csv("/Users/eliasshaeffer/Desktop/somatic_gatk_low/summaryWithVariants.csv")