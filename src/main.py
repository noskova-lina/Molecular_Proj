import os
from database import create_database, insert_xyz_data, view_table, delete_database

db_path = 'xyz_molecules.db'
# Absolute path to the directory with the current script (src)
current_directory = os.path.dirname(os.path.abspath(__file__))
# Path to the project root directory
project_directory = os.path.dirname(current_directory)
# Relative path to the data folder
directory_path = os.path.join(project_directory, 'data', 'structures')

create_database(db_path)
insert_xyz_data(db_path, directory_path)
view_table(db_path)
# delete_database(db_path)
