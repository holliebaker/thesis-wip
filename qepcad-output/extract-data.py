#!/usr/bin/python3

import sys
import re
from tabulate import tabulate

r = 0

# regexes
label_regex = re.compile("^([A-Z])_([0-9]+),([0-9]+)\s+=\s((?:input  \*\*\* Refinement of cell (\([0-9,]\)) \*\*\*  )|(.*))$")
# groups in order: type, level, index, parent 
polynomial_regex = re.compile("^\s+=\s(.*)$")
# one group, matching polynomial

cell_start_regex = re.compile("^-* Information.*\(([0-9,]+)\)")
# contains the cell index
cell_end_regex = re.compile("^-+$")
# line of dashes, indicates end of cell
dimension_regex = re.compile("^Dimension\s*:\s*([0-9]+)")
# dimension
signpf_regex = re.compile("^Level ([0-9]+)\s*:\s\(([+0\-,]*)\)")
# level %d : list of (+,-,0)
alpha_regex = re.compile("(?:^alpha = .*|^)the unique root of ([0-9+\-^x\s]+)\sbetween ([0-9\-\/]+) and ([0-9\-\/]+)")
# polynomial, lower and upper
coordinate_regex = re.compile("^Coordinate ([0-9]+)\s*=\s*(.*)$")
# k, value
mathmode_regex = re.compile("\d+,\d+|fac|res|dis|input|alpha")
# matches index, so we can wrap in curlies
def mathmode_regex_fun(matches):
    # extract text
    text = matches.group(0)
    is_op = True

    # special case for alpha
    if text == "alpha":
        return "\\" + text

    # wrap in curlies
    text = "{" + text + "}"

    if text[1].isdigit():
        # it's a comma-separated list of numbers
        return text

    # otherwise it's akeyword
    return "\mathrm" + text


# constants
PROJECTION_FACTORS = "d-proj-fac"
PROJECTION_POLYNOMIALS = "d-proj-pol"
INPUT_POLYNOMIALS = "d-input"
REFINEMENT_POLYNOMIALS = "d-ref"
TRUE_CELLS = "d-true-cells"

# data
projection_factors = {}
projection_polynomials = {}
input_polynomials = {}
refinement_polynomials = {}
cells = []

def process_alpha(matches):
    # only include alpha if it is a polynomial
    a = matches.group(2)
    b = matches.group(3)
    if a == b:
        return None

    alpha = {
        "minimal": matches.group(1),
        "interval": (matches.group(2), matches.group(3))
    }
     
    return alpha

# process polynomials written in qepcad output
def process_polynomials(lines):
    global r
    data = {}
    level = 0
    label = None
    parent = ""
    polynomial = None
    alpha = None
    point = {}
    
    empty_lines = 0
    while lines != []:
        line = lines.pop(0)
        # track number of empty lines
        if line == "":
            empty_lines += 1

            if empty_lines > 0 and level > 0:
                if level not in data.keys():
                    r = max(level,r)
                    data[level] = []
                
                data[level].append({
                    "label": label,
                    "parent": parent,
                    "polynomial": polynomial,
                    "alpha": alpha,
                    "point": point
                })

                level = 0
                label = None
                parent = ""
                polynomial = None
                alpha = None
                point = {}

            if empty_lines == 3:
                # 3 empty lines indicates end of output
                return data
            else:
                continue

        # otherwise, line not empty
        empty_lines = 0
        
        # do a regex match
        matches = label_regex.search(line)
        if matches:
            # parse label and parent
            level = int(matches.group(2))
            label = {
                "type": matches.group(1),
                "level": level,
                "index": int(matches.group(3))
            }

            if matches.group(5):
                parent = matches.group(5)
            else:
                parent = matches.group(4)
                        
        matches = polynomial_regex.search(line)
        if matches:
            polynomial = matches.group(1)
            
            if point or polynomial.startswith("The"):
                # it's a sample point
                polynomial = None

            continue

        matches = alpha_regex.search(line)
        if matches:
            alpha = process_alpha(matches)

            continue

        # point type, parse coordinate            
        matches = coordinate_regex.search(line)
        if matches:
            k = int(matches.group(1))
            point[k] = matches.group(2)

            continue
        
    return data

def map_signpf(s):
    if s == "-":
        return -1
    elif s == "+":
        return 1
    else:
        return 0

# process cad cells
def process_cells(lines):
    data = []
    index = None
    dimension = 0
    signpf = {}
    alpha = None
    sample = {}
    
    while lines != []:
        line = lines.pop(0)
        
        # end of all cells
        if line.startswith("Before"):
            return data

        # end of one cell
        matches = cell_end_regex.search(line)
        if matches and index:
            data.append({
                "index": index,
                "dimension": dimension,
                "signpf": signpf,
                "alpha": alpha,
                "sample": sample
            })
            
            index = None
            dimension = 0
            signpf = {}
            alpha = None
            sample = {}
            continue

        # start of one cell, find index
        matches = cell_start_regex.search(line)
        if matches:
            index = list(map(int, matches.group(1).split(",")))
            continue
            
        # sign of projection factor, level k
        matches = signpf_regex.search(line)
        if matches:
            k = int(matches.group(1))
            signpf[k] = list(map(map_signpf, matches.group(2).split(",")))
            continue

        matches = dimension_regex.search(line)
        if matches:
            dimension = int(matches.group(1))
            continue

        matches = alpha_regex.search(line)
        if matches:
            # only include alpha if it is a polynomial
            a = matches.group(2)
            b = matches.group(3)
            if not a == b:
                alpha = {
                    "minimal": matches.group(1),
                    "interval": (matches.group(2), matches.group(3))
                }
            
            continue

        matches = coordinate_regex.search(line)
        if matches:
            k = int(matches.group(1))
            sample[k] = matches.group(2)

            if sample[k].startswith("the"): # extended sample point, same format as alpha
                sample[k] = process_alpha(alpha_regex.search(sample[k]))
                
                if not sample[k]:
                    sample[k] = "0"

            continue

    return data            

# print polynomials
       
# open and read qepcad file
filename = sys.argv[1]
with open(filename) as f:
    lines = f.read().splitlines()

# read lines
while lines != []:
    line = lines.pop(0)
    if line.startswith(PROJECTION_FACTORS):
        projection_factors = process_polynomials(lines)
    elif line.startswith(PROJECTION_POLYNOMIALS):
        projection_polynomials = process_polynomials(lines)
    elif line.startswith(INPUT_POLYNOMIALS):
        input_polynomials = process_polynomials(lines)
    elif line.startswith(REFINEMENT_POLYNOMIALS):
        print("refinement polys")
        refinement_polynomials = process_polynomials(lines)
    elif line.startswith(TRUE_CELLS):
        cells = process_cells(lines)

def mathmode(text):
    text = mathmode_regex.sub(mathmode_regex_fun, text)
    return "\[" + text + "\]"


def polynomial_format_row(poly, ptype):
    label = poly["label"]
    label_text = mathmode(label["type"] + "_" + str(label["level"]) + "," + str((label["index"])))
    alpha = poly["alpha"]
    point = poly["point"]
    parent = poly["parent"]

    if label["type"] == "Q":
        ptype = ptype + ", **QA**"
        
    if ptype == "RP" and label["type"] == "M":
        ptype = ptype + ", " + mathmode(parent)
        parent = None
    
    if parent and not parent.startswith("input"):
        parent = mathmode(parent)
    else:
        parent = ""

    polynomial = poly["polynomial"]
    if not polynomial:
        # sample point type
        polynomial = mathmode(
            ",\> ".join([mathmode(point[k]) for k in point])
        )

        if alpha:
            parent = mathmode(alpha["minimal"] + ",\>(" + cell["alpha"]["interval"][0] + ", " + cell["alpha"]["interval"][1] + ")")
    else:
        # standard
        polynomial = mathmode(polynomial)

    return {
        "Label": label_text,
        "Type": ptype,
        "Polynomial": polynomial,
        "Parent": parent,
    }


def flatten_polynomials(poly_map):
    global r
    flattened = []

    # level
    for k in range(1,r+1):
        # type of poly
        for t in poly_map:
            polys = poly_map[t]
            
            # type has level
            if k in polys:
                flattened.extend(
                    [polynomial_format_row(x, t) for x in polys[k]]
                )
            
    return flattened

sign_operator_map = {
    -1: " < ",
    0: " = ",
    1: " > "
}

def cell_format_row(cell):
    global projection_factors
    alpha = cell["alpha"]
    if alpha:
        alpha = mathmode(alpha["minimal"] + ",\>(" + cell["alpha"]["interval"][0] + ", " + cell["alpha"]["interval"][1] + ")")
    else:
        alpha = "Rational"

    # first row: 
    rows = [{
        "Position": mathmode(", ".join(list(map(str, cell["index"])))),
        "Dimension": mathmode(str(cell["dimension"])),
        "k": "",
        "Sample": alpha,
        "Signs of projection factors": ""
    }]
    
    # sample and signs of projection factors
    for k in range(1,r+1):
        sign_proj_factors = []
        pfs = projection_factors[k]
        signs = cell["signpf"][k]

        for poly in pfs:
            if poly["label"]["type"] == "M":
                # refinement polynomial
                continue
                
            i = poly["label"]["index"] - 1
            sign_proj_factors.append(poly["polynomial"] + sign_operator_map[signs[i]] + "0")
    
        rows.append({
            "Position": "",
            "Dimension": "",
            "k": mathmode(str(k)),
            "Sample": mathmode(cell["sample"][k]),
            "Signs of projection factors": mathmode(",\> ".join(sign_proj_factors))
        })
    
    return rows

polys_table = tabulate(
    flatten_polynomials({
        "IP": input_polynomials,
        "PP": projection_polynomials,
        "PF": projection_factors,
        "RP": refinement_polynomials
    }),
    headers="keys",
    tablefmt="github"
)

cells_table = tabulate(
    [row for cell in cells for row in cell_format_row(cell)],
    headers="keys",
    tablefmt="github"
)

print("**Polynomials**")
print(polys_table)
print()
print("**Cells**")
print(cells_table)


