import numpy as np
import sys
import os

def read_bestsqs(file):
    try:
        with open(file, 'r') as f:
            lines = f.readlines()
        
        dArrVec1 = np.array([list(map(float, lines[i].split())) for i in range(3)])
        dArrVec2 = np.array([list(map(float, lines[i+3].split())) for i in range(3)])
        
        dArrAtom, strArrAtom = [], []
        for line in lines[6:]:
            x, y, z, elem = line.split()
            dArrAtom.append([float(x), float(y), float(z)])
            strArrAtom.append(elem)
        
        nAtom = len(dArrAtom)
        dArrAtom = np.array(dArrAtom)
        
        # Count number of elements
        unique_elements = sorted(set(strArrAtom), key=strArrAtom.index)
        nElem = len(unique_elements)
        nArrElem = [strArrAtom.count(elem) for elem in unique_elements]
        return True, dArrVec1, dArrVec2, dArrAtom, strArrAtom, nAtom, nElem, nArrElem, unique_elements
    except Exception as e:
        print(f"Error reading file {file}: {e}")
        return False, None, None, None, None, None, None, None, None

def inverse_matrix(m1):
    try:
        inv = np.linalg.inv(m1)
        return True, inv
    except np.linalg.LinAlgError:
        return False, None

def main(file):
    success, dArrVec1, dArrVec2, dArrAtom, strArrAtom, nAtom, nElem, nArrElem, strArrElem = read_bestsqs(file)
    if not success:
        print(f"Errors occurred when reading {file}!")
        sys.exit(1)
    
    dArrLatVec = dArrVec2 @ dArrVec1
    success, dArrLatVecInv = inverse_matrix(dArrLatVec)
    if not success:
        print("Errors occurred when converting to direct.")
        sys.exit(1)
    
    dArrAtom2 = dArrAtom @ dArrLatVecInv.T
    
    output_file = f"{file}-POSCAR"
    with open(output_file, 'w') as f:
        f.write("POSCAR\n")
        f.write("1.0\n")
        for row in dArrLatVec:
            f.write(" ".join(f"{val:12.8f}" for val in row) + "\n")
        f.write(" ".join(strArrElem) + "\n")
        f.write(" ".join(map(str, nArrElem)) + "\n")
        f.write("Direct\n")
        for row in dArrAtom2:
            f.write(" ".join(f"{val:12.8f}" for val in row) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: sqs2poscar.py FileName")
        sys.exit(1)
    main(sys.argv[1])
