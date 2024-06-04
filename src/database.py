import sqlite3
import os
import glob
import logging
from settings import db_path, directory_path

def create_database():
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    # Create the molecules table if it does not exist
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

def insert_xyz_data(directory_path):
    file_paths = sorted(glob.glob(os.path.join(directory_path, '*.xyz')))
    insert_xyz_files(file_paths)

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

    conn.close()
    return data
