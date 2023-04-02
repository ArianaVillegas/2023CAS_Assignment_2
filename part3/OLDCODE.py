#OLD CODE
'''
syn_count = 0
nonsyn_count = 0
temp = []
temp2 = []
temp3 = []
checked = []

for amino, acid in base_muts.mutations.items():
    if amino in checked:
        continue
    syn_count += len(base_muts.mutations[amino][0])
    nonsyn_count += len(base_muts.mutations[amino][1])
    for amino2 in base_muts.mutations.keys():
        if amino2 == amino:
            continue
        if len(acid[0]) != 0:
            if len(acid[0]) == 3:
                if acid[0][0] == amino2 or acid[0][1] == amino2 or acid[0][2] == amino2:
                    syn_count += 3
                    nonsyn_count += len(acid[1])
                    checked.append(amino2)
            elif len(acid[0]) == 2:
                if acid[0][0] == amino2 or acid[0][1] == amino2:
                    syn_count += 2
                    nonsyn_count += len(acid[1])
                    checked.append(amino2)
            elif len(acid[0]) == 1:
                if acid[0] == amino2:
                    syn_count += 1
                    nonsyn_count += len(acid[1])
                    checked.append(amino2)
        else:
            continue
    print(syn_count)
    temp.append(syn_count)
    temp2.append(nonsyn_count+syn_count)
    temp3.append([syn_count, nonsyn_count+syn_count])
    if len(temp) != 0:
        if syn_count != temp[-1]:
            temp.append(syn_count)
    else:
        temp.append(syn_count)

    #syn_count = 0
    #nonsyn_count = 0

for amino in table.values():
    for amino2 in table.values():
        if amino == amino2:
            syn_count += 1
        else:
            nonsyn_count += 1
    print(syn_count, nonsyn_count)
    syn_count = 0
    nonsyn_count = 0
'''
