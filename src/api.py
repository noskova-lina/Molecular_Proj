from fastapi import FastAPI, HTTPException
import os

from database import create_database, insert_xyz_data, view_table, delete_database

app = FastAPI()

db_path = 'xyz_molecules.db'
directory_path = './data/structures'

@app.get("/")
def read_root():
    return {"Hello": "Welcome to the Molecular Data API"}

@app.post("/create-database/")
def create_db():
    try:
        create_database(db_path)
        return {"message": "Database created successfully"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/insert-data/")
def insert_data():
    try:
        insert_xyz_data(db_path, directory_path)
        return {"message": "Data inserted successfully"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/view-data/")
def get_data():
    try:
        data = view_table(db_path)
        return data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.delete("/delete-database/")
def delete_db():
    try:
        delete_database(db_path)
        return {"message": "Database deleted successfully"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
