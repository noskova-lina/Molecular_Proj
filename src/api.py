from fastapi import FastAPI, HTTPException, Query, Body
from settings import db_path, directory_path
from database import insert_xyz_file, insert_xyz_files, insert_xyz_data, view_table

app = FastAPI()


@app.get("/")
def read_root():
    return {"Hello": "Welcome to the Molecular Data API"}


@app.post("/insert-data/")
def insert_data():
    try:
        insert_xyz_data(directory_path)
        return {"message": "Data inserted successfully"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/insert-file/")
def insert_file(xyz_file_path: str):
    try:
        insert_xyz_file(xyz_file_path)
        return {"message": f"Data from {xyz_file_path} inserted successfully"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/insert-files/")
def insert_files(xyz_file_paths: list[str] = Body(...)):
    try:
        insert_xyz_files(xyz_file_paths)
        return {"message": "Data from files inserted successfully"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/view-data/")
def get_data(molecule_name: str = Query(default=None)):
    try:
        if molecule_name is None:
            raise HTTPException(status_code=400, detail="Molecule name must be provided")

        data = view_table(molecule_name)
        if not data:
            raise HTTPException(status_code=404, detail=f"Molecule '{molecule_name}' not found")

        return data
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
