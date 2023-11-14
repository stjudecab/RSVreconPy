import os
import re

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
                    cov = a + c + g + t

                    if a / cov > 0.5:
                        sequence += "A"
                    elif c / cov > 0.5:
                        sequence += "C"
                    elif g / cov > 0.5:
                        sequence += "G"
                    elif t / cov > 0.5:
                        sequence += "T"
                    else:
                        if a > c and a > g and a > t:
                            sequence += "A"
                        elif c > a and c > g and c > t:
                            sequence += "C"
                        elif g > a and g > c and g > t:
                            sequence += "G"
                        elif t > a and t > g and t > c:
                            sequence += "T"
                else:
                    print("error")
            else:
                print("error1")

    return sequence

# determine the subtype
def determine_subtype(subtype_dict):
    sorted_subtypes = sorted(subtype_dict.keys(), key=lambda k: subtype_dict[k], reverse=True)
    return sorted_subtypes[0]
