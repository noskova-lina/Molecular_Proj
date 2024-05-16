from database import create_database, insert_xyz_data, view_table, delete_database

db_path = 'xyz_molecules.db'
directory_path = 'C:/Users/AlinaNoskova/PycharmProjects/Molecular_Proj/data/structures'

create_database(db_path)
insert_xyz_data(db_path, directory_path)
view_table(db_path)
#delete_database(db_path)
