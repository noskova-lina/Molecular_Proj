import uvicorn
from settings import directory_path
from database import create_database, insert_xyz_data, insert_xyz_files, view_table, delete_database
import logging


create_database()
insert_xyz_data(directory_path)
# view_table('dsgdb9nsd_000002')
# delete_database()

if __name__ == "__main__":
    uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=True)


def main():
    logging.basicConfig(level=logging.INFO)

    try:
        create_database()
        logging.info("Database created successfully.")
    except Exception as e:
        logging.error(f"Failed to create database: {str(e)}")

    file_paths = [
        '../data/sample.xyz',
        '../data/sample_2.xyz'
    ]

    try:
        insert_xyz_files(file_paths)
        logging.info("Files inserted into the database.")
    except Exception as e:
        logging.error(f"Failed to insert files: {str(e)}")

    try:
        data = view_table('sample')
        if data:
            logging.info("Data retrieved successfully.")
            for row in data:
                print(row)
        else:
            logging.error("Molecule not found.")
    except Exception as e:
        logging.error(f"Failed to retrieve data: {str(e)}")

    try:
        delete_database()
        logging.info("Database deleted successfully.")
    except Exception as e:
        logging.error(f"Failed to delete database: {str(e)}")

# if __name__ == "__main__":
#     main()
