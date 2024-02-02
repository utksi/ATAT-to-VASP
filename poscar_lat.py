import numpy as np
import sys

def read_poscar(poscar_file):
    with open(poscar_file, 'r') as file:
        lines = file.readlines()

    scaling_factor = float(lines[1].strip())
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]]) * scaling_factor
    atom_types = lines[5].strip().split()
    atom_counts = list(map(int, lines[6].strip().split()))
    coordinate_type = lines[7].strip()
    positions_start_line = 8 if coordinate_type.lower().startswith(('c', 'd')) else 7
    atomic_positions = np.array([list(map(float, line.split()[:3])) for line in lines[positions_start_line:positions_start_line+sum(atom_counts)]])

    return lattice_vectors, atom_types, atom_counts, atomic_positions, coordinate_type.lower().startswith('d')

def lattice_parameters(lattice_vectors):
    a = np.linalg.norm(lattice_vectors[0])
    b = np.linalg.norm(lattice_vectors[1])
    c = np.linalg.norm(lattice_vectors[2])
    alpha = np.arccos(np.dot(lattice_vectors[1], lattice_vectors[2]) / (b * c)) * 180 / np.pi
    beta = np.arccos(np.dot(lattice_vectors[0], lattice_vectors[2]) / (a * c)) * 180 / np.pi
    gamma = np.arccos(np.dot(lattice_vectors[0], lattice_vectors[1]) / (a * b)) * 180 / np.pi
    return a, b, c, alpha, beta, gamma

def write_lat_in(filename, lattice_vectors, atom_types, atom_counts, atomic_positions):
    a, b, c, alpha, beta, gamma = lattice_parameters(lattice_vectors)
    with open(filename, 'w') as file:
        file.write(f"{a} {b} {c} {alpha} {beta} {gamma}  # Coordinate system: a b c alpha beta gamma\n")
        file.write("1 0 0  # Primitive unit cell\n0 1 0\n0 0 1\n")
        for atom_type, count in zip(atom_types, atom_counts):
            file.write(f"0 0 0 {atom_type}  # Atoms in the lattice\n")
            if count > 1:
                for i in range(1, count):
                    file.write(f"{' '.join(map(str, atomic_positions[i-1]))} {atom_type}\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python poscar_lat.py POSCAR [lat.in]")
        sys.exit(1)

    poscar_file = sys.argv[1]
    lat_in_file = 'lat.in' if len(sys.argv) == 2 else sys.argv[2]

    lattice_vectors, atom_types, atom_counts, atomic_positions, is_direct = read_poscar(poscar_file)
    if not is_direct:
        print("Warning: Atomic positions are not in direct coordinates. Manual conversion may be required.")
    write_lat_in(lat_in_file, lattice_vectors, atom_types, atom_counts, atomic_positions)
    print(f"Converted {poscar_file} to {lat_in_file} in ATAT format.")

if __name__ == "__main__":
    main()
