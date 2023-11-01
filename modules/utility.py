def select_subset(data,x,y,stride=1):
    #Select subset of simulations in plaquette data
    keys = [key for key in data]
    keys = keys[x:y:stride]
    subset_data = {}
    try:
        for key in keys:
            subset_data[key] = data[key]
    except:
        subset_data[keys] = data[keys]
        
    return subset_data