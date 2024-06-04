import requests
import pytest

BASE_URL = "http://localhost:8000"

def test_read_root():
    response = requests.get(f"{BASE_URL}/")
    assert response.status_code == 200
    assert response.json() == {"Hello": "Welcome to the Molecular Data API"}

def test_insert_data():
    response = requests.post(f"{BASE_URL}/insert-data/")
    assert response.status_code == 200
    assert "Data inserted successfully" in response.json()["message"]

def test_insert_file():
    xyz_file_path = "../data/sample.xyz"
    response = requests.post(f"{BASE_URL}/insert-file/", data={"xyz_file_path": xyz_file_path})
    assert response.status_code == 200
    assert f"Data from {xyz_file_path} inserted successfully" in response.json()["message"]

def test_get_data():
    response = requests.get(f"{BASE_URL}/view-data/")
    assert response.status_code == 200
    assert isinstance(response.json(), list)