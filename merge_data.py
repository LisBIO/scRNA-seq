import pandas as pd

def mysplit(x):
    a = x.split(".")[0]
    return a


import os
data = pd.read_csv("./ERR1211/count_file/ERR1211178_Count.txt", sep = "\t", header=None)
data = data[:55487]
gene_name = data[0].apply(mysplit)
# print(data[0])


filelist = os.listdir("./ERR1211/count_file/")
new = pd.DataFrame(index=data.index)
# print(new)
filelist.remove(".DS_Store")
filelist.remove(".txt")
# print(len(filelist))
for file in filelist:
    print(file)
    data = pd.read_csv("./ERR1211/count_file/" + file, sep="\t", header = None)
    str1 = file.split(".")
    # print(str1)
    cell = str1[0].split("_")[0]
    # if len(str1) > 3:
    #     str2 = str1[0].split("_")
    #     if len(str2) == 2:
    #         cell = str1[0].split("_")[1] + "_" + str1[1]
    #     else:
    #         cell = str2[2] + "_" + str2[1] + "_" + str1[1]
    # else:
    #     cell = str1[0].split("_")[1] + "_" + str1[0].split("_")[2]

    data[cell] = data[1]
    new = new.merge(data[[cell]], right_index=True, left_index=True)
    print("Finish" + " " + cell)

new.index = gene_name
new.to_csv("./mouse.csv")
""""""
