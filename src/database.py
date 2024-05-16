import sqlite3
import os
import glob


def create_database(db_path):
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

    conn.close()


def insert_xyz_data(db_path, directory_path):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    for xyz_file_path in sorted(glob.glob(os.path.join(directory_path, '*.xyz'))):
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


def view_table(db_path):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute('SELECT * FROM molecules ORDER BY molecule_name, atom_index LIMIT 19')

    column_names = [description[0] for description in cur.description]
    print(column_names)

    rows = cur.fetchall()
    for row in rows:
        print(row)

    conn.close()


def delete_database(db_path):
    try:
        os.remove(db_path)
        print("Database deleted successfully.")
    except OSError as e:
        print(f"Error: {e.strerror}. File {e.filename} could not be deleted.")
