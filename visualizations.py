import matplotlib.pyplot as plt

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import confusion_matrix
import sortedcontainers



def create_heatmap_from_data(data, title="Comparison between trajectories / clusterings", names=None, heatmap_label='Trajectory Nr.'):
    import numpy as np
    import matplotlib.pyplot as plt

    num_individuals = max(max(row[0] for row in data), max(row[1] for row in data)) + 1

    if names is None:
        names = [str(i) for i in range(num_individuals)]

    heatmap_data = np.zeros((num_individuals, num_individuals))

    for row in data:
        i, j, value = row
        heatmap_data[i][j] = value

    plt.figure(figsize=(8, 6))
    plt.imshow(heatmap_data, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Values')
    plt.xlabel(heatmap_label)
    plt.ylabel(heatmap_label)
    plt.title(title)
    
    plt.xticks(ticks=np.arange(num_individuals), labels=names, rotation=90)
    plt.yticks(ticks=np.arange(num_individuals), labels=names)

    plt.show()



def get_trajectory_heatmap(clusterings, clusterings2, nmi = True, create_heatmap = True, title = "Comparison between trajectories / clusterings", names=None, heatmap_label='Trajectory Nr.'):
    import sklearn
    rand_array = []
    for i1, cl1 in  enumerate(clusterings):
        for i2, cl2 in enumerate(clusterings2):

            if not nmi:
                rand_array.append([i1,i2, sklearn.metrics.adjusted_rand_score(cl1, cl2)])
            else:
                rand_array.append([i1,i2, sklearn.metrics.normalized_mutual_info_score(cl1, cl2)])

    
    if create_heatmap:
        if nmi:
            create_heatmap_from_data(rand_array, title= title + "\n NMI", heatmap_label=heatmap_label, names=names)
        else:
            create_heatmap_from_data(rand_array, title= title + "\n adj. Rand", heatmap_label=heatmap_label, names=names)


    return rand_array    




def plot_residue_line(residues, labels):
    labels_set = sortedcontainers.SortedSet(labels)
    num_labels = len(labels_set)
    
    color_map = plt.cm.get_cmap('tab10', num_labels) 

    fig, ax = plt.subplots()

    for i in range(len(residues) - 1):
        x = [residues[i], residues[i + 1]]
        y = [0, 0]  
        label_color = color_map(labels_set.index(labels[i]))
        ax.plot(x, y, color=label_color, linewidth=14)  

    plt.axis('off')

    ax.annotate(f'{residues[0]}', (residues[0], 0), xytext=(5, -15), textcoords='offset points', ha='center', color='black')
    ax.annotate(f'{residues[-1]}', (residues[-1], 0), xytext=(-5, -15), textcoords='offset points', ha='center', color='black')


    plt.text((residues[0] + residues[-1]) / 2, -0.01, 'Residue number', ha='center')

    plt.show()



def plot_residue_line_plus_noise(residues, labels, noise_label = -1, title = None):
    labels_set = sortedcontainers.SortedSet(labels)
    num_labels = len(labels_set)
    
    color_map = plt.cm.get_cmap('tab10', num_labels) 

    fig, ax = plt.subplots()

    for i in range(len(residues) - 1):
        x = [residues[i], residues[i + 1]]
        y = [0, 0]  
        label_color = color_map(labels_set.index(labels[i]))
        if noise_label:
            if labels[i] == noise_label:
                label_color = "black"
        ax.plot(x, y, color=label_color, linewidth=14)  

    plt.axis('off')

    ax.annotate(f'{residues[0]}', (residues[0], 0), xytext=(5, -15), textcoords='offset points', ha='center', color='black')
    ax.annotate(f'{residues[-1]}', (residues[-1], 0), xytext=(-5, -15), textcoords='offset points', ha='center', color='black')


    plt.text((residues[0] + residues[-1]) / 2, -0.01, 'Residue number', ha='center')

    if title:
        plt.title(title)

    plt.show()


def multiple_line_plots(data, items=None, titles=None, full_title="Residue clustering", orig_data = None, noise_value = -1):

    import sortedcontainers
    labels_set = [y for y in [(item['labels']) for item in data] for y in y] 
    #print(labels_set)
    labels_set = sortedcontainers.SortedSet(labels_set)
    #print(labels_set)


    num_labels = len(labels_set)
    color_map = plt.cm.get_cmap('tab20', num_labels)  

    num_plots = len(data)
    num_cols = 2  
    num_rows = (num_plots + num_cols - 1) // num_cols  


    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 2 * num_rows), 
                            gridspec_kw={'hspace': 0.00001})  # Adjust the value as needed

    add_counter = 0

    if num_rows > 1:
        for i, ax_row in enumerate(axs):

            
            for j, ax in enumerate(ax_row):
                index = i * num_cols + j

                if index < num_plots:

                    item_data = data[index]

                    items = item_data['items']
                    labels = item_data['labels']
                    
                    for k in range(len(items) - 1):
                        x = [items[k], items[k + 1]]
                        y = [0, 0]  
                        label_color = color_map(labels_set.index(labels[k]))

                        if orig_data != None:
                            if orig_data[add_counter][k] == noise_value:
                                label_color = 'black'
     
                        ax.plot(x, y, color=label_color, linewidth=24)  
                        ax.axis('off') 

                    
                    ax.annotate(f'{items[0]}', (items[0], 0), xytext=(5, -26), textcoords='offset points', ha='center', color='black')
                    ax.annotate(f'{items[-1]}', (items[-1], 0), xytext=(-5, -26), textcoords='offset points', ha='center', color='black')

                    
                    ax.text((items[0] + items[-1]) / 2, -0.03, 'Residue number', ha='center')

                    ax.set_title(item_data['title'])

                else:
                    ax.remove()
                
                add_counter = add_counter + 1
    else:


        index = 0

        while index < num_plots :
            print("index")
            print(index)
            if index < num_plots :
                ax = axs[index]

                item_data = data[index]

                items = item_data['items']
                labels = item_data['labels']
                
                for k in range(len(items) - 1):
                    x = [items[k], items[k + 1]]
                    y = [0, 0] 
                    label_color = color_map(labels_set.index(labels[k]))
                    
                    ax.plot(x, y, color=label_color, linewidth=24)  
                    ax.axis('off')  

                ax.annotate(f'{items[0]}', (items[0], 0), xytext=(5, -26), textcoords='offset points', ha='center', color='black')
                ax.annotate(f'{items[-1]}', (items[-1], 0), xytext=(-5, -26), textcoords='offset points', ha='center', color='black')

                ax.text((items[0] + items[-1]) / 2, -0.03, 'Residue number', ha='center')

                ax.set_title(item_data['title'])
            index = index + 1

    plt.suptitle(full_title)

    plt.tight_layout()

    plt.show()





def rearrange_labels_multi(*label_lists):
    """
    Rearrange multiple lists of labels to match the first list using the Hungarian (Munkres) algorithm.
    """
    num_lists = len(label_lists)
    num_labels = max(max(labels) for labels in label_lists) + 1
    
    agg_cm = np.zeros((num_labels, num_labels))
    
    for labels1 in label_lists:
        for labels2 in label_lists:
            cm = confusion_matrix(labels1, labels2, labels=range(num_labels))
            agg_cm += cm
    
    row_ind, col_ind = linear_sum_assignment(-agg_cm)  
    
    rearranged_label_lists = [[col_ind[label] for label in labels] for labels in label_lists]
    
    return rearranged_label_lists



def create_list_of_dicts(labels, items=None, titles=None):
    if items is None:
        items = []
        for label in labels:
            print(label)
            items.append(list(range(len(label))))

    if titles is None:
        titles = []
        for il, label in enumerate(labels):
            titles.append("plot nr. " + str(il))


    list_of_dicts = []

    for i in range(len(labels)):
        title = titles[i]
        item = items[i]
        label = labels[i]
        try:
            data_dict = {'title': title, 'items': item, 'labels': label.tolist()}
        except:
            data_dict = {'title': title, 'items': item, 'labels': label}
        list_of_dicts.append(data_dict)
    return list_of_dicts


def reassign_labels(cluster_labels1, cluster_labels2):
    unique_labels1 = np.unique(cluster_labels1)
    unique_labels2 = np.unique(cluster_labels2)
    
    num_clusters1 = len(unique_labels1)
    num_clusters2 = len(unique_labels2)
    
    cost_matrix = np.zeros((num_clusters1, num_clusters2))
    for i in range(num_clusters1):
        for j in range(num_clusters2):
            common_elements = np.sum((cluster_labels1 == unique_labels1[i]) & (cluster_labels2 == unique_labels2[j]))
            cost_matrix[i, j] = -common_elements 
    
    # Use the Hungarian (Munkres) algorithm to find the optimal assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    label_mapping = {}
    for i, j in zip(row_ind, col_ind):
        if cost_matrix[i, j] != -np.inf: 
            label_mapping[unique_labels2[j]] = unique_labels1[i]
    
    reassigned_labels2 = [label_mapping[label] if label in label_mapping else label for label in cluster_labels2]
    
    return cluster_labels1, reassigned_labels2



def iterate_rearrange_labels(labels_lists):
    print(labels_lists)
    try:
     start_label = labels_lists[0].tolist()
    except:
     start_label = labels_lists[0]
       
    #print(start_label)
    rearranged_lists = []
    rearranged_lists.append(start_label)

    for cluster_list in labels_lists[1:]:

          if isinstance(cluster_list, list):
               cluster_list = cluster_list[0]

          new_cluster_list = reassign_labels(start_label, cluster_list)
          print(new_cluster_list)
          rearranged_lists.append(new_cluster_list[-1])

    return rearranged_lists



def line_plot_workflow(data, titles = None, full_title = "residue line plots", rearrange = True, hdb_scan_noise=False):
    from copy import deepcopy
    data_final = deepcopy(data)
    if rearrange:
        data_final = iterate_rearrange_labels(data)
    dicts = create_list_of_dicts(data_final, titles = titles)

    if not hdb_scan_noise:
        multiple_line_plots(dicts, full_title = full_title)
    else:
        multiple_line_plots(dicts, full_title = full_title, orig_data = data)


def find_noise_resiudes(clusters, noise_nr=-1):
    indices = [i for i, x in enumerate(clusters) if x == noise_nr]
    return indices

