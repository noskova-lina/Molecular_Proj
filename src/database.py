import csv
import sqlite3
import os
import logging
import zipfile

from rdkit import Chem
from rdkit.Geometry import Point3D
from settings import db_path


def create_database():
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute('''
    CREATE TABLE IF NOT EXISTS molecules (
        molecule_name TEXT NOT NULL,
        atom_index INTEGER NOT NULL,
        atom TEXT NOT NULL,
        x REAL NOT NULL,
        y REAL NOT NULL,
        z REAL NOT NULL,
        PRIMARY KEY (molecule_name, atom_index)
    )''')

    cur.execute('CREATE INDEX IF NOT EXISTS idx_molecule_name ON molecules (molecule_name)')

    cur.execute('''
    CREATE TABLE IF NOT EXISTS scalar_coupling_constants (
        id INTEGER PRIMARY KEY,
        molecule_name TEXT NOT NULL,
        atom_index_0 INTEGER NOT NULL,
        atom_index_1 INTEGER NOT NULL,
        type TEXT NOT NULL,
        scalar_coupling_constant REAL NOT NULL,
        FOREIGN KEY (molecule_name) REFERENCES molecules(molecule_name)
    )''')

    cur.execute('''
    CREATE TABLE IF NOT EXISTS scalar_coupling_contributions (
        molecule_name TEXT NOT NULL,
        atom_index_0 INTEGER NOT NULL,
        atom_index_1 INTEGER NOT NULL,
        type TEXT NOT NULL,
        fc REAL NOT NULL,
        sd REAL NOT NULL,
        pso REAL NOT NULL,
        dso REAL NOT NULL,
        FOREIGN KEY (molecule_name) REFERENCES molecules(molecule_name)
    )''')

    conn.close()


def delete_database():
    try:
        os.remove(db_path)
    except OSError as e:
        logging.error(f"Error: {e.strerror}. File {e.filename} could not be deleted.")


def insert_xyz_file(xyz_file_path):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    with open(xyz_file_path, 'r') as file:
        lines = file.readlines()
        molecule_name = os.path.basename(xyz_file_path).replace('.xyz', '')
        for index, line in enumerate(lines[2:]):
            parts = line.strip().split()
            atom_type = parts[0]
            x, y, z = map(float, parts[1:])
            cur.execute(
                'INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z) VALUES (?, ?, ?, ?, ?, ?)',
                (molecule_name, index, atom_type, x, y, z))

    conn.commit()
    conn.close()


def insert_xyz_files(file_paths):
    for xyz_file_path in file_paths:
        insert_xyz_file(xyz_file_path)


def insert_xyz_data(zip_file_path):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        # List all the files in the zip archive
        xyz_files = [f for f in zip_ref.namelist() if f.endswith('.xyz')]

        for xyz_file_name in sorted(xyz_files):
            with zip_ref.open(xyz_file_name) as file:
                lines = file.readlines()
                molecule_name = os.path.basename(xyz_file_name).replace('.xyz', '')
                for index, line in enumerate(lines[2:]):
                    parts = line.decode('utf-8').strip().split()
                    atom_type = parts[0]
                    x, y, z = map(float, parts[1:])
                    cur.execute(
                        'INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z) VALUES (?, ?, ?, ?, ?, ?)',
                        (molecule_name, index, atom_type, x, y, z))

    conn.commit()
    conn.close()


def view_table(molecule_name=None):
    if not molecule_name:
        logging.error("Molecule name is not provided.")
        return []

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute('SELECT * FROM molecules WHERE molecule_name = ? ORDER BY molecule_name, atom_index', (molecule_name,))
    rows = cur.fetchall()

    if not rows:
        logging.error(f"Molecule '{molecule_name}' not found in the database.")
        return []

    column_names = [description[0] for description in cur.description]
    data = [dict(zip(column_names, row)) for row in rows]
    print(data)

    conn.close()
    return data


def get_molecule_data(molecule_name):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute('''
    SELECT * FROM molecules WHERE molecule_name = ?
    ''', (molecule_name,))
    rows = cur.fetchall()
    conn.close()
    return rows


def insert_scalar_coupling_constants_from_zip(zip_file_scalar_coupling_path: str):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    with zipfile.ZipFile(zip_file_scalar_coupling_path, 'r') as zip_ref:
        # List all the files in the zip archive
        csv_files = [f for f in zip_ref.namelist() if f.endswith('.csv') and 'scalar_coupling_constants' in f]

        for csv_file_name in csv_files:
            with zip_ref.open(csv_file_name) as csvfile:
                reader = csv.DictReader(csvfile.read().decode('utf-8').splitlines())
                for row in reader:
                    cur.execute('''
                    INSERT INTO scalar_coupling_constants (id, molecule_name, atom_index_0, atom_index_1, type, scalar_coupling_constant)
                    VALUES (?, ?, ?, ?, ?, ?)
                    ''', (
                        row['id'],
                        row['molecule_name'],
                        row['atom_index_0'],
                        row['atom_index_1'],
                        row['type'],
                        row['scalar_coupling_constant']
                    ))

    conn.commit()
    conn.close()


def insert_scalar_coupling_contributions_from_zip(zip_file_scalar_contributions_path: str):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    with zipfile.ZipFile(zip_file_scalar_contributions_path, 'r') as zip_ref:
        # List all the files in the zip archive
        csv_files = [f for f in zip_ref.namelist() if f.endswith('.csv') and 'scalar_coupling_contributions' in f]

        for csv_file_name in csv_files:
            with zip_ref.open(csv_file_name) as csvfile:
                reader = csv.DictReader(csvfile.read().decode('utf-8').splitlines())
                for row in reader:
                    cur.execute('''
                    INSERT INTO scalar_coupling_contributions (molecule_name, atom_index_0, atom_index_1, type, fc, sd, pso, dso)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (
                        row['molecule_name'],
                        row['atom_index_0'],
                        row['atom_index_1'],
                        row['type'],
                        row['fc'],
                        row['sd'],
                        row['pso'],
                        row['dso']
                    ))

    conn.commit()
    conn.close()


def get_molecule_data(molecule_name):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute('''
    SELECT * FROM molecules WHERE molecule_name = ?
    ''', (molecule_name,))
    rows = cur.fetchall()
    conn.close()
    return rows


def get_scalar_coupling_constants(molecule_name):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute('''
    SELECT * FROM scalar_coupling_constants WHERE molecule_name = ?
    ''', (molecule_name,))
    rows = cur.fetchall()
    conn.close()
    return rows


def get_scalar_coupling_contributions(molecule_name):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute('''
    SELECT * FROM scalar_coupling_contributions WHERE molecule_name = ?
    ''', (molecule_name,))
    rows = cur.fetchall()
    conn.close()
    return rows


def get_data_by_atom_index_and_type(molecule_name: str, atom_index_0: int, atom_index_1: int, type: str):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute('''
    SELECT * FROM scalar_coupling_constants
    WHERE molecule_name = ? AND atom_index_0 = ? AND atom_index_1 = ? AND type = ?
    ''', (molecule_name, atom_index_0, atom_index_1, type))
    rows = cur.fetchall()
    conn.close()
    return rows


def get_molecules_by_atom(atom: str):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute('''
    SELECT DISTINCT molecule_name FROM molecules WHERE atom = ?
    ''', (atom,))
    rows = cur.fetchall()
    conn.close()
    return [row[0] for row in rows]


def get_molecules_by_coordinate_range(x_min: float, x_max: float, y_min: float, y_max: float, z_min: float,
                                      z_max: float):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute('''
    SELECT DISTINCT molecule_name FROM molecules
    WHERE x BETWEEN ? AND ? AND y BETWEEN ? AND ? AND z BETWEEN ? AND ?
    ''', (x_min, x_max, y_min, y_max, z_min, z_max))
    rows = cur.fetchall()
    conn.close()
    return [row[0] for row in rows]


def get_xyz_data(molecule_name: str):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM molecules WHERE molecule_name = ?", (molecule_name,))
    data = cursor.fetchone()
    conn.close()
    return data


def update_molecule_data(molecule_name: str, atom_index: int, atom: str = None, x: float = None, y: float = None,
                         z: float = None):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    columns = []
    values = []
    if atom is not None:
        columns.append("atom = ?")
        values.append(atom)
    if x is not None:
        columns.append("x = ?")
        values.append(x)
    if y is not None:
        columns.append("y = ?")
        values.append(y)
    if z is not None:
        columns.append("z = ?")
        values.append(z)
    values.append(molecule_name)
    values.append(atom_index)

    query = f'''
    UPDATE molecules
    SET {", ".join(columns)}
    WHERE molecule_name = ? AND atom_index = ?
    '''
    cur.execute(query, tuple(values))
    conn.commit()
    conn.close()


def xyz_to_rdkit_molecule(xyz_data: str):
    lines = xyz_data.strip().split("\n")
    num_atoms = int(lines[0])
    mol = Chem.RWMol()

    for line in lines[2:2 + num_atoms]:
        parts = line.split()
        atom_symbol = parts[0]
        x, y, z = map(float, parts[1:4])
        atom = Chem.Atom(atom_symbol)
        idx = mol.AddAtom(atom)
        mol.GetConformer().SetAtomPosition(idx, Point3D(x, y, z))

    Chem.AllChem.EmbedMolecule(mol)
    return mol
