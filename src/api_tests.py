import requests
import time

BASE_URL = "http://localhost:8000"


def test_read_root():
    start_time = time.time()
    response = requests.get(f"{BASE_URL}/")
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json() == {"Hello": "Welcome to the Molecular Data API"}
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_insert_data():
    start_time = time.time()
    response = requests.post(f'{BASE_URL}/insert-data/')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json().get("message") == "Data inserted successfully"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_insert_file():
    start_time = time.time()
    file_path = '../data/structure_sample.xyz'
    response = requests.post(f'{BASE_URL}/insert-file/?xyz_file_path={file_path}')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json().get("message") == f"Data from {file_path} inserted successfully"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_insert_files():
    start_time = time.time()
    file_paths = [
        '../data/structure_sample_2.xyz',
        '../data/structure_sample_3.xyz'
    ]
    response = requests.post(f'{BASE_URL}/insert-files/', json=file_paths)
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json().get("message") == "Data from files inserted successfully"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_insert_scalar_coupling_constants_file():
    start_time = time.time()
    file_path = '../data/scalar_coupling_constants.zip'
    response = requests.post(f'{BASE_URL}/insert-scalar-coupling-constants-file/', json={'file_path': file_path})
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json().get("message") == "Scalar coupling constants inserted successfully from file"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_insert_scalar_coupling_contributions_file():
    start_time = time.time()
    file_path = '../data/scalar_coupling_contributions.zip'
    response = requests.post(f'{BASE_URL}/insert-scalar-coupling-contributions-file/', json={'file_path': file_path})
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json().get("message") == "Scalar coupling contributions inserted successfully from file"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_view_data():
    start_time = time.time()
    molecule_name = 'structure_sample'
    response = requests.get(f'{BASE_URL}/view-data/?molecule_name={molecule_name}')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json(), "Expected non-empty response"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_view_data_invalid_name():
    start_time = time.time()
    invalid_name = 'non_existent_molecule'
    response = requests.get(f'{BASE_URL}/view-data/?molecule_name={invalid_name}')
    assert response.status_code == 404, f"Expected status code 404 but got {response.status_code}. Response: {response.text}"
    assert response.json().get("detail") == f"Molecule '{invalid_name}' not found"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_get_molecule_data():
    start_time = time.time()
    molecule_name = 'dsgdb9nsd_000001'
    response = requests.get(f'{BASE_URL}/get-molecule-data/?molecule_name={molecule_name}')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json(), "Expected non-empty response"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_get_scalar_coupling_constants():
    start_time = time.time()
    molecule_name = 'dsgdb9nsd_000001'
    response = requests.get(f'{BASE_URL}/get-scalar-coupling-constants/?molecule_name={molecule_name}')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json(), "Expected non-empty response"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_get_scalar_coupling_contributions():
    start_time = time.time()
    molecule_name = 'dsgdb9nsd_000001'
    response = requests.get(f'{BASE_URL}/get-scalar-coupling-contributions/?molecule_name={molecule_name}')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json(), "Expected non-empty response"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_get_data_by_atom_index_and_type():
    start_time = time.time()
    molecule_name = 'dsgdb9nsd_000001'
    atom_index_0 = 1
    atom_index_1 = 0
    type = '1JHC'
    response = requests.get(
        f'{BASE_URL}/get-data-by-atom-index-and-type/?molecule_name={molecule_name}&atom_index_0={atom_index_0}&atom_index_1={atom_index_1}&type={type}')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json(), "Expected non-empty response"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_get_molecules_by_atom():
    start_time = time.time()
    atom = 'H'
    response = requests.get(f'{BASE_URL}/get-molecules-by-atom/?atom={atom}')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json(), "Expected non-empty response"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_get_molecules_by_coordinate_range():
    start_time = time.time()
    x_min, x_max = 0.0, 1.0
    y_min, y_max = 0.0, 1.0
    z_min, z_max = 0.0, 1.0
    response = requests.get(
        f'{BASE_URL}/get-molecules-by-coordinate-range/?x_min={x_min}&x_max={x_max}&y_min={y_min}&y_max={y_max}&z_min={z_min}&z_max={z_max}')
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json(), "Expected non-empty response"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")


def test_update_molecule_data():
    start_time = time.time()
    update_data = {
        'molecule_name': 'dsgdb9nsd_000001',
        'atom_index': 0,
        'x': 1.234
    }
    response = requests.put(
        f'{BASE_URL}/update-molecule-data/',
        json=update_data
    )
    assert response.status_code == 200, f"Expected status code 200 but got {response.status_code}. Response: {response.text}"
    assert response.json().get("message") == "Molecule data updated successfully"
    duration = time.time() - start_time
    print(f"Test duration: {duration:.4f} seconds")
