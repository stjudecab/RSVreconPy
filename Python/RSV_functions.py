import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

# get sub folders for a dir
def get_sub_folders(folder_path):
    sub_folders = []

    with os.scandir(folder_path) as entries:
        for entry in entries:
            if entry.is_dir() and not entry.name.startswith(('.', '..')):
                sub_folders.append(entry.name)

    return sub_folders

# test if elements in array A in array B, if not, return elements that are not in B
def elements_not_in_array(array_a, array_b):
    not_in_b = [element for element in array_a if element not in array_b]
    return not_in_b

# sum the pct values
def pct_sum(*args):
    total = sum(float(item.replace('%', '')) for item in args)
    return total

# assemble sequences from IGV counts
def processIGV(file_name):
    sequence = ""

    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('track') or line.startswith('#') or line.startswith('variableStep chrom='):
                continue
            elif line[0].isdigit():
                match = re.match(r'^(\d+)\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0.+', line)
                if match:
                    pos, a, c, g, t = map(int, match.groups())
                    
                    count_dict = {'A':a,'C':c,'G':g,'T':t}
                    cov = a + c + g + t

                    if cov > 0:
                        max_key = max(count_dict, key=count_dict.get)
                        sequence += max_key
                    else:
                        sequence += "N"
                else:
                    print("error")
            else:
                print("error1")

    return sequence

# determine the subtype
def determine_subtype(subtype_dict):
    sorted_subtypes = sorted(subtype_dict.keys(), key=lambda k: subtype_dict[k], reverse=True)
    if subtype_dict[sorted_subtypes[0]] > 30:
        return sorted_subtypes[0]
    elif subtype_dict[sorted_subtypes[0]] > 10:
        return sorted_subtypes[0].'(Low mapping rate)'
    else:
        return 'Not RSV'

# check tool availability
def check_tool_availability(tool_name):
    try:
        # Use subprocess to run the command with the "--version" flag
        subprocess.run([tool_name, '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"{tool_name} is available.")
        return 0
    except:
        print(f"{tool_name} is not available. Please install it and try again.")
        return 1

# Python code to get difference of two lists not using set()
def list_diff(li1, li2):
    li_dif = [i for i in li1 if i not in li2]
    return li_dif

# generate coverage fig from wig file
def generate_Cov_fig(wig_file, fig_file):
    # read wig file and transfer it to a data frame
    data = {}
    last_pos = 1

    with open(wig_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('track') or line.startswith('#') or line.startswith('variableStep chrom='):
                continue
            elif line[0].isdigit():
                match = re.match(r'^(\d+)\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0.+', line)
                if match:
                    pos, a, c, g, t = map(int, match.groups())
                    data[pos] = [a, c, g, t]
                    last_pos = pos

    index_all = list(range(1,last_pos))
    cur_index = list(data.keys())
    empty_index = list_diff(index_all, cur_index)

    if len(empty_index) > 0:
        for index in empty_index:
            data[index] = [0,0,0,0]

    df = pd.DataFrame(data, index=['A', 'C', 'G', 'T'])

    # generate coverage plot
    ax = df.transpose().plot(kind='bar', stacked=True, colormap='viridis')
    ax.set_xlabel('Position', fontsize=18)
    ax.set_ylabel('Count (Log scale)', fontsize=18)
    ax.set_title('Stacked Bar Plot of A, T, C, G Counts at Different Positions')
    ax.xaxis.label.set_size(20)
    plt.yscale('symlog')

    # save plot
    plt.gcf().set_size_inches(30, 3)
    plt.savefig(fig_file, dpi = 300)

    # return results
    return df