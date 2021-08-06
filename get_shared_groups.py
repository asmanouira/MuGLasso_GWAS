def merge_LD_groups(groups_list_pop1,groups_list_pop2):
    shared_groups = [1]
    for i in range(1,len(groups_list_pop1)):

        if groups_list_pop1[i] == groups_list_pop2[i] == 1:
            shared_groups.append(groups_list_pop1[i])
        elif groups_list_pop1[i] == groups_list_pop2[i] and groups_list_pop1[i-1] == groups_list_pop1[i] and groups_list_pop2[i-1] == groups_list_pop2[i]:
            shared_groups.append(shared_groups[i-1])    
        elif groups_list_pop1[i] != groups_list_pop2[i] and groups_list_pop1[i-1] == groups_list_pop1[i] and groups_list_pop2[i-1] == groups_list_pop2[i]:
            shared_groups.append(shared_groups[i-1])
        elif groups_list_pop1[i] == groups_list_pop2[i] <= shared_groups[i-1]:
            shared_groups.append(shared_groups[i-1]+1)
        elif groups_list_pop1[i] != groups_list_pop2[i] and (groups_list_pop1[i-1] != groups_list_pop1[i] or groups_list_pop2[i-1]!= groups_list_pop2[i]):
            shared_groups.append( shared_groups[i-1]+1)
    return(shared_groups)