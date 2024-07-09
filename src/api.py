from fastapi import FastAPI, HTTPException, Query, Body
from settings import directory_path
from database import insert_xyz_file, insert_xyz_files, insert_xyz_data, view_table, get_xyz_data, xyz_to_rdkit_molecule
from pydantic import BaseModel
from typing import Optional
from io import BytesIO
from fastapi.responses import StreamingResponse
from rdkit import Chem

from src.database import (
    get_scalar_coupling_contributions,
    get_scalar_coupling_constants,
    get_molecule_data,
    insert_scalar_coupling_contributions_from_csv,
    insert_scalar_coupling_constants_from_csv,
    get_data_by_atom_index_and_type,
    get_molecules_by_atom,
    get_molecules_by_coordinate_range,
    update_molecule_data
)

app = FastAPI()


class MoleculeUpdate(BaseModel):
    molecule_name: str
    atom_index: int
    atom: Optional[str] = None
    x: Optional[float] = None
    y: Optional[float] = None
    z: Optional[float] = None


@app.get("/")
def read_root():
    return {"Hello": "Welcome to the Molecular Data API"}


@app.post("/insert-data/")
def insert_data():
    try:
        insert_xyz_data(directory_path)
        return {"message": "Data inserted successfully"}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Directory not found")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/insert-file/")
def insert_file(xyz_file_path: str):
    try:
        insert_xyz_file(xyz_file_path)
        return {"message": f"Data from {xyz_file_path} inserted successfully"}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"File {xyz_file_path} not found")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/insert-files/")
def insert_files(xyz_file_paths: list[str] = Body(...)):
    try:
        insert_xyz_files(xyz_file_paths)
        return {"message": "Data from files inserted successfully"}
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/insert-scalar-coupling-constants-file/")
def insert_scalar_coupling_constants_file(file_path: str = Body(..., embed=True)):
    try:
        insert_scalar_coupling_constants_from_csv(file_path)
        return {"message": "Scalar coupling constants inserted successfully from file"}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"File {file_path} not found")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/insert-scalar-coupling-contributions-file/")
def insert_scalar_coupling_contributions_file(file_path: str = Body(..., embed=True)):
    try:
        insert_scalar_coupling_contributions_from_csv(file_path)
        return {"message": "Scalar coupling contributions inserted successfully from file"}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"File {file_path} not found")
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


@app.get("/get-molecule-data/")
def get_molecule_data_api(molecule_name: str):
    try:
        data = get_molecule_data(molecule_name)
        if not data:
            raise HTTPException(status_code=404, detail=f"No data found for molecule '{molecule_name}'")
        return data
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/get-scalar-coupling-constants/")
def get_scalar_coupling_constants_api(molecule_name: str):
    try:
        data = get_scalar_coupling_constants(molecule_name)
        if not data:
            raise HTTPException(status_code=404,
                                detail=f"No scalar coupling constants found for molecule '{molecule_name}'")
        return data
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/get-scalar-coupling-contributions/")
def get_scalar_coupling_contributions_api(molecule_name: str):
    try:
        data = get_scalar_coupling_contributions(molecule_name)
        if not data:
            raise HTTPException(status_code=404,
                                detail=f"No scalar coupling contributions found for molecule '{molecule_name}'")
        return data
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/get-data-by-atom-index-and-type/")
def get_data_by_atom_index_and_type_api(molecule_name: str, atom_index_0: int, atom_index_1: int, type: str):
    try:
        data = get_data_by_atom_index_and_type(molecule_name, atom_index_0, atom_index_1, type)
        if not data:
            raise HTTPException(status_code=404,
                                detail=f"No data found for molecule '{molecule_name}' with atom_index_0 {atom_index_0}, atom_index_1 {atom_index_1}, and type {type}")
        return data
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/get-molecules-by-atom/")
def get_molecules_by_atom_api(atom: str):
    try:
        data = get_molecules_by_atom(atom)
        if not data:
            raise HTTPException(status_code=404, detail=f"No molecules found containing atom '{atom}'")
        return data
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/get-molecules-by-coordinate-range/")
def get_molecules_by_coordinate_range_api(x_min: float, x_max: float, y_min: float, y_max: float, z_min: float,
                                          z_max: float):
    try:
        data = get_molecules_by_coordinate_range(x_min, x_max, y_min, y_max, z_min, z_max)
        if not data:
            raise HTTPException(status_code=404, detail=f"No molecules found within the coordinate range")
        return data
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.put("/update-molecule-data/")
def update_molecule_data_api(update: MoleculeUpdate):
    try:
        update_molecule_data(update.molecule_name, update.atom_index, update.atom, update.x, update.y, update.z)
        return {"message": "Molecule data updated successfully"}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/visualize-molecule-from-db/")
def visualize_molecule_from_db(molecule_name: str):
    try:
        xyz_data = get_xyz_data(molecule_name)
        if not xyz_data:
            raise HTTPException(status_code=404, detail="Molecule not found in the database")

        mol = xyz_to_rdkit_molecule(xyz_data)
        if mol is None:
            raise HTTPException(status_code=400, detail="Unable to create molecule from XYZ data")

        img = Chem.Draw.MolToImage(mol)
        img_byte_arr = BytesIO()
        img.save(img_byte_arr, format='PNG')
        img_byte_arr.seek(0)

        return StreamingResponse(img_byte_arr, media_type="image/png")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
