network = open("summary_frac_tp_100_perm_1--frac_fp_0_perm_1_1.tsv", "r") #input file
pred_file = open("pred_groups.txt", "r") #file that has pred.group.[num] and respective TF names
TF_network = open("TF_network.txt", "w") #output file

pred = pred_file.readlines()
list_pred = []
for i in range(len(pred)):
    line = pred[i]
    line = line.split()
    list_pred.append(line)

header = network.readline()
TF_network.write(header)
data = network.readlines()
final = []

for i in range(len(data)):
    line = data[i]
    line = line.split()
    if line[0] == 'pred.group.3':
        for j in range(len(list_pred[2])-1):
            L = []
            TF = list_pred[2][j+1]
            L.append(TF)
            s = "\t".join(line[1:])
            L.append(s)
            s1 = '\t'.join(L)
            # print(s1)
            TF_network.write(s1)
            TF_network.write("\n")
            
    elif line[0] == 'pred.group.2':
        for j in range(len(list_pred[1])-1):
            L = []
            TF = list_pred[1][j+1]
            L.append(TF)
            s = "\t".join(line[1:])
            L.append(s)
            s1 = '\t'.join(L)
            # print(s1)
            TF_network.write(s1)
            TF_network.write("\n")
    
    elif line[0] == 'pred.group.1':
        for j in range(len(list_pred[0])-1):
            L = []
            TF = list_pred[0][j+1]
            L.append(TF)
            s = "\t".join(line[1:])
            L.append(s)
            s1 = '\t'.join(L)
            # print(s1)
            TF_network.write(s1)
            TF_network.write("\n")
    elif line[0] == 'pred.group.4':
        for j in range(len(list_pred[3])-1):
            L = []
            TF = list_pred[3][j+1]
            L.append(TF)
            s = "\t".join(line[1:])
            L.append(s)
            s1 = '\t'.join(L)
            # print(s1)
            TF_network.write(s1)
            TF_network.write("\n")
    elif line[0] == 'pred.group.5':
        for j in range(len(list_pred[4])-1):
            L = []
            TF = list_pred[4][j+1]
            L.append(TF)
            s = "\t".join(line[1:])
            L.append(s)
            s1 = '\t'.join(L)
            # print(s1)
            TF_network.write(s1)
            TF_network.write("\n")
    elif line[0] == 'pred.group.6':
        for j in range(len(list_pred[5])-1):
            L = []
            TF = list_pred[5][j+1]
            L.append(TF)
            s = "\t".join(line[1:])
            L.append(s)
            s1 = '\t'.join(L)
            # print(s1)
            TF_network.write(s1)
            TF_network.write("\n")
    elif line[0] == 'pred.group.7':
        for j in range(len(list_pred[6])-1):
            L = []
            TF = list_pred[6][j+1]
            L.append(TF)
            s = "\t".join(line[1:])
            L.append(s)
            s1 = '\t'.join(L)
            # print(s1)
            TF_network.write(s1)
            TF_network.write("\n")
#     elif line[0] == 'pred.group.8':
#         for j in range(len(list_pred[7])-1):
#             L = []
#             TF = list_pred[7][j+1]
#             L.append(TF)
#             s = "\t".join(line[1:])
#             L.append(s)
#             s1 = '\t'.join(L)
#             # print(s1)
#             TF_network.write(s1)
#             TF_network.write("\n")
#     elif line[0] == 'pred.group.9':
#         for j in range(len(list_pred[8])-1):
#             L = []
#             TF = list_pred[8][j+1]
#             L.append(TF)
#             s = "\t".join(line[1:])
#             L.append(s)
#             s1 = '\t'.join(L)
#             # print(s1)
#             TF_network.write(s1)
#             TF_network.write("\n")
    else:
        s = "\t".join(line[0:])
        TF_network.write(s)
        TF_network.write("\n")